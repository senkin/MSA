#!/usr/bin/env python3
import argparse
import json
import os
import sys
import warnings
import pandas as pd
import numpy as np
import time
from scipy.optimize import nnls
from numba import jit
from common_methods import make_folder_if_not_exists, calculate_similarity, format_threshold

# Vectorized bootstrap methods
BOOTSTRAP_METHODS = ['binomial', 'multinomial', 'multinomial_weight', 'poisson', 'classic',
                     'bootstrap_residuals']


def resample_channels(n_channels):
    """Row positions of a classic (i.i.d., with replacement) resample of n_channels rows."""
    # same draw as DataFrame.sample(n, replace=True) on the global RNG
    return np.random.choice(n_channels, size=n_channels, replace=True)


def bootstrap_mutation_table(input_dataframe, method="classic", fitted=None, residuals=None, rows=None):
    """One bootstrap resample of a (channels x samples) mutation table.

    classic: rows resampled with replacement (i.i.d. assumption); the resampled rows keep
        the original channel labels, as rows are exchangeable under this assumption.
    bootstrap_residuals: semi-parametric: fitted + residuals resampled across channels as
        in classic, negative counts set to 0. input_dataframe is not used.
    rows: row positions to use for classic / bootstrap_residuals (default: a new draw), so
        several tables can share one resample.
    """
    if method in ("classic", "bootstrap_residuals"):
        if method == "bootstrap_residuals":
            if residuals is None or fitted is None:
                raise ValueError(f"Residuals and fitted values required for {method}")
            table = residuals
        else:
            table = input_dataframe
        if rows is None:
            rows = resample_channels(len(table))
        resampled = pd.DataFrame(table.values[rows], index=table.index, columns=table.columns)
        if method == "classic":
            return resampled
        return (fitted + resampled.values).clip(lower=0)

    elif method == "poisson":
        # Vectorized Poisson sampling
        return input_dataframe * np.random.poisson(1, size=input_dataframe.shape)

    elif method in ["binomial", "multinomial", "multinomial_weight"]:
        sums = input_dataframe.sum(axis=0).values
        norm_mutations = input_dataframe.div(sums, axis=1).fillna(0).values
        bootstrap_data = np.zeros_like(input_dataframe.values)

        for i, col_sum in enumerate(sums):
            if method == "binomial":
                bootstrap_data[:, i] = np.random.binomial(int(col_sum), norm_mutations[:, i])
            elif method == "multinomial":
                bootstrap_data[:, i] = np.random.multinomial(int(col_sum), norm_mutations[:, i])
            elif method == "multinomial_weight":
                n = len(norm_mutations[:, i])
                counts = np.random.multinomial(n, [1/n] * n)
                bootstrap = np.repeat(input_dataframe.values[:, i], counts)
                np.random.shuffle(bootstrap)
                bootstrap_data[:, i] = bootstrap

        return pd.DataFrame(bootstrap_data, index=input_dataframe.index, columns=input_dataframe.columns)

    else:
        raise ValueError(f"Unknown bootstrap method: {method}")


@jit(nopython=True, cache=True)
def calculate_stats_numba(observed, fitted):
    """JIT-compiled statistics calculation for speed."""
    residuals = observed - fitted
    rss = np.sum(residuals ** 2)
    chi2 = np.sum(residuals ** 2 / (fitted + 1e-10))  # Avoid division by zero

    mean_obs = np.mean(observed)
    tss = np.sum((observed - mean_obs) ** 2)
    r2 = 1 - rss / tss if tss > 0 else 0

    return rss, chi2, r2


# Mapping from a stat_info index (see the layout assembled in
# perform_signature_attribution) to the (metric, normalise) pair understood by
# calculate_similarity. Index 8 (L2_normalised_by_first) is handled directly for
# speed and is the metric used by the optimisation loops throughout the pipeline.
_IDX_TO_METRIC = {
    4: ('Cosine', False),
    5: ('Correlation', False),
    6: ('Chebyshev', True),
    7: ('L1', True),
    9: ('L3', True),
    10: ('jensen-shannon', False),
}


def nnls_batched(A, B, masks, b_index=None, out_weights=None, out_fitted=None):
    """Solve, for each problem ``j``, NNLS(A[:, masks[:, j]], B[:, b_index[j]]).

    This is the general batched primitive and the CPU implementation behind
    batched_solve_and_score(). It is NOT the GPU swap point: the optimisation hot
    path never consumes ``out_weights`` or ``out_fitted``, only the single
    similarity scalar that batched_solve_and_score() reduces them to. Returning
    the full fitted matrix costs n_channels * 8 bytes per problem to deliver 8
    bytes of signal (an n_channels:1 amplification), so a GPU backend should fuse
    the reduction and implement batched_solve_and_score() instead.

    This function is still needed by the callers that genuinely want weights:
    batched_nnls_shared_A() (non-optimised attribution).

    Parameters
    ----------
    A : ndarray (n_channels, n_signatures)
        Full signature matrix (shared by all problems).
    B : ndarray (n_channels, n_targets)
        Distinct target vectors. Targets are *not* duplicated per problem; the
        ``b_index`` gather keeps memory at O(n_targets) rather than O(n_problems).
    masks : ndarray (n_signatures, n_problems), bool
        Column j selects which columns of A are active for problem j.
    b_index : ndarray (n_problems,), int, optional
        Target column of B for each problem. Defaults to the identity mapping
        (which requires n_problems == n_targets).
    out_weights : ndarray (n_signatures, n_problems), optional
        Dense weights, zero outside the mask.
    out_fitted : ndarray (n_channels, n_problems), optional

    Returns
    -------
    (out_weights, out_fitted)
    """
    n_problems = masks.shape[1]
    if b_index is None:
        b_index = np.arange(n_problems)
    if out_weights is None:
        out_weights = np.zeros((A.shape[1], n_problems))
    if out_fitted is None:
        out_fitted = np.zeros((A.shape[0], n_problems))

    for j in range(n_problems):
        cols = np.flatnonzero(masks[:, j])
        w, _ = nnls(A[:, cols], B[:, b_index[j]])
        out_weights[:, j] = 0.0
        out_weights[cols, j] = w
        out_fitted[:, j] = A[:, cols] @ w
    return out_weights, out_fitted


def batched_nnls_shared_A(A, B):
    """Solve NNLS(A, b) for every column ``b`` of ``B`` (shared design matrix).

    Used by the non-optimised attribution path, where every sample is fitted
    against the same signature matrix with no column masking. This is the
    all-columns-active special case of nnls_batched() and delegates to it.

    Unlike the optimisation hot path, this caller genuinely needs the weights, so
    it cannot go through batched_solve_and_score(). It is also a minor workload -
    one run per dataset and mutation type - so it is a low-value GPU target.

    Should it ever be worth accelerating: with a constant ``A`` and no masking, a
    factorisation of ``A`` can be computed once and reused across all N solves -
    impossible in the masked case, where every problem selects a different subset
    of columns. Specialise this function if so, not nnls_batched().

    Parameters
    ----------
    A : ndarray (n_channels, n_signatures)
    B : ndarray (n_channels, n_samples)

    Returns
    -------
    ndarray (n_signatures, n_samples) : weights for each sample (column).
    """
    masks = np.ones((A.shape[1], B.shape[1]), dtype=bool)
    weights, _ = nnls_batched(A, B, masks)
    return weights


def _l2_similarity_batched(targets, fitted, norm_targets):
    """Vectorised L2_normalised_by_first similarity for a batch of fits.

    Batched counterpart of _similarity_from_fitted() for idx == 8:
        1 - ||target - fitted|| / ||target||
    """
    return 1.0 - np.linalg.norm(targets - fitted, axis=0) / norm_targets


# GPU backend for batched_solve_and_score(). None keeps the pure-CPU/scipy path;
# main() replaces it with a _GpuBatchedSolver when --use_gpu is passed. Kept as a
# module global so the deep greedy call sites need no extra plumbing.
_GPU_SOLVER = None

# Process-global GPU environment, set up once by _init_gpu_backend(). The stream is
# held here so it is not garbage-collected while it is the active cupy stream.
_GPU_INITIALISED = False
_GPU_STREAM = None
# Device-memory ceiling for GPU allocations (the RMM pool cap); suggest_chunk() sizes
# the NNLS batch against this.
_GPU_POOL_BYTES = 0

# Auto chunk-size knobs (used only when --gpu_batch_size <= 0). _CHUNK_MEM_FRACTION is
# the share of the RMM pool one batched call's transient allocations may use (the rest
# covers the resident A/B/||b|| copies, cuML overhead and fragmentation);
# _CHUNK_WORKSPACE_FACTOR is the headroom multiplier on the estimate of cuML's
# per-problem Lawson-Hanson workspace, which we cannot measure directly. Raise the
# fraction / lower the factor to push bigger batches once calibrated on a real card.
_CHUNK_MEM_FRACTION = 0.5
_CHUNK_WORKSPACE_FACTOR = 3.0
_CHUNK_MIN = 256
_CHUNK_MAX = 1 << 20


def _init_gpu_backend():
    """Bring up an RMM pool and put cupy on raft's stream (once per process).

    Mirrors the reference GPU setup from the cuML batched-NNLS PR (rapidsai/cuml#8402):

    * A bounded RMM pool routes cupy's allocations so the greedy's many short-lived
      device buffers do not thrash cudaMalloc/cudaFree. The cap comes from the
      launcher via RMM_POOL_SIZE_MIB (MiB, already divided across the processes
      sharing the GPU); absent that it defaults to 60% of total memory. It grows on
      demand, not up front.
    * cupy is switched onto the per-thread default stream (ptds), the same stream
      raft's Handle() uses. Sharing one stream keeps cupy's H2D copies / gathers
      ordered before the cuML kernels that read them; on separate streams that race
      surfaces intermittently as CUDA_ERROR_ILLEGAL_ADDRESS at handle.sync().
    """
    global _GPU_INITIALISED, _GPU_STREAM, _GPU_POOL_BYTES
    if _GPU_INITIALISED:
        return
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator
    import cupy as cp

    env_mib = os.environ.get("RMM_POOL_SIZE_MIB")
    if env_mib is not None:
        pool_size = int(float(env_mib) * (1 << 20))
    else:
        _, total_memory = cp.cuda.runtime.memGetInfo()
        pool_size = int(total_memory * 0.6)
    pool_size -= pool_size % 256  # RMM requires the pool size to be a multiple of 256 bytes
    _GPU_POOL_BYTES = pool_size
    rmm.reinitialize(pool_allocator=True, initial_pool_size=0, maximum_pool_size=pool_size)
    cp.cuda.set_allocator(rmm_cupy_allocator)

    _GPU_STREAM = cp.cuda.Stream(ptds=True)
    _GPU_STREAM.use()
    _GPU_INITIALISED = True


class _GpuBatchedSolver:
    """cuML-backed implementation of batched_solve_and_score().

    Wraps ``cuml.solvers.nnls_batched`` (rapidsai/cuml PR #8402), whose signature
    and semantics match our own nnls_batched(): ``nnls_batched(A, B, masks,
    b_index=..., compute_fitted=True) -> (X, fitted)``, with A shared, B gathered by
    b_index and masks column-major (n_signatures, n_problems).

    Follows the GPU setup from the cuML batched-NNLS PR (rapidsai/cuml#8402): a bounded RMM pool and cupy
    on raft's per-thread default stream (see _init_gpu_backend), plus A / B / ||b||
    kept device-resident for as long as the same (A, B, norm_obs) objects keep
    arriving - they are loop-invariant across the greedy steps of one attribution,
    so each uploads at most once per attribution. The residual-norm reduction is
    fused on the device, so only the (n_problems,) similarity vector returns to the
    host (fitted never leaves the GPU).

    Precision (--gpu_precision, see configure_gpu_backend):
      'fp64'  - solve and score in fp64 (matches the example run on our float64
                arrays; closest to the CPU reference).
      'fp32'  - solve and score in fp32 (fastest; may flip a small fraction of
                borderline selections vs the CPU reference).
      'mixed' - solve in fp32 (the throughput win) but score in fp64: the greedy
                compares tiny similarity differences (base_sim - trial_sim) that are
                cancellation-prone, so B / ||b|| are also kept in fp64 for the
                reduction. Default.
    """

    def __init__(self, precision='mixed'):
        _init_gpu_backend()
        from cuml.solvers import nnls_batched as cuml_nnls_batched
        self._nnls_batched = cuml_nnls_batched

        if precision == 'fp64':
            self._solve_dtype, self._score_dtype = np.float64, np.float64
        elif precision == 'fp32':
            self._solve_dtype, self._score_dtype = np.float32, np.float32
        elif precision == 'mixed':
            self._solve_dtype, self._score_dtype = np.float32, np.float64
        else:
            raise ValueError(f"Unknown GPU precision '{precision}' (use mixed, fp32 or fp64)")
        self._mixed = self._solve_dtype != self._score_dtype
        self._cache = {}
        self.precision = precision
        self.auto_batch_size = None

    def provenance(self):
        """What produced a GPU result: device, driver/runtime, library versions, the
        cuML build (commit recorded by setup_gpu_env.sh), precision and batch sizing.

        Everything is best-effort: a missing field is reported, never raised.
        """
        info = {'precision': self.precision,
                'rmm_pool_bytes': _GPU_POOL_BYTES,
                'auto_batch_size': self.auto_batch_size}
        try:
            import cupy as cp
            import rmm
            import cuml
            info['versions'] = {'cuml': getattr(cuml, '__version__', None),
                                'cupy': getattr(cp, '__version__', None),
                                'rmm': getattr(rmm, '__version__', None)}
            dev_id = cp.cuda.Device().id
            props = cp.cuda.runtime.getDeviceProperties(dev_id)
            name = props.get('name', b'')
            info['device'] = {'id': dev_id,
                              'name': name.decode() if isinstance(name, bytes) else str(name),
                              'compute_capability': f"{props.get('major')}.{props.get('minor')}",
                              'total_memory_bytes': props.get('totalGlobalMem')}
            info['cuda'] = {'driver_version': cp.cuda.runtime.driverGetVersion(),
                            'runtime_version': cp.cuda.runtime.runtimeGetVersion()}
        except Exception as e:  # provenance must never break a run
            info['error'] = f'{type(e).__name__}: {e}'
        build_file = os.path.join(sys.prefix, 'etc', 'msa_cuml_build.json')
        try:
            with open(build_file) as f:
                info['cuml_build'] = json.load(f)
        except OSError:
            info['cuml_build'] = f'not recorded ({build_file} missing; see setup_gpu_env.sh)'
        return info

    def _resident(self, A, B, norm_obs):
        """Device copies of the loop-invariant A, B, ||b||, cached by identity.

        A new (A, B, norm_obs) triple evicts the previous generation so device
        memory stays bounded; the host objects are pinned in the cache so their
        id() cannot be reused while cached. In 'mixed' mode B is kept in both the
        solve dtype (fp32) and the score dtype (fp64).
        """
        import cupy as cp
        gen = (id(A), id(B), id(norm_obs))
        if self._cache.get("gen") != gen:
            self._cache.clear()
            self._cache["gen"] = gen
            self._cache["refs"] = (A, B, norm_obs)
            self._cache["A"] = cp.asfortranarray(cp.asarray(A, dtype=self._solve_dtype))
            self._cache["B_solve"] = cp.asfortranarray(cp.asarray(B, dtype=self._solve_dtype))
            self._cache["B_score"] = (self._cache["B_solve"] if not self._mixed
                                      else cp.asarray(B, dtype=self._score_dtype))
            self._cache["norm"] = cp.asarray(norm_obs, dtype=self._score_dtype)
        return (self._cache["A"], self._cache["B_solve"],
                self._cache["B_score"], self._cache["norm"])

    def suggest_chunk(self, n_channels, n_signatures):
        """Auto-size the NNLS batch (problems per cuML call) to the GPU's memory.

        Picks the largest chunk whose estimated device footprint fits a fraction of the
        RMM pool, so a bigger card / smaller context (e.g. SBS-96) gets a larger batch
        and a memory-heavy context (SBS-1536/4608) a smaller one, with no manual tuning.
        Used only when --gpu_batch_size <= 0. See the _CHUNK_* knobs above.
        """
        solve_bytes = np.dtype(self._solve_dtype).itemsize
        score_bytes = np.dtype(self._score_dtype).itemsize
        # Buffers we allocate per problem (exact): fitted + gathered target + residual
        # over n_channels; masks (uint8) + discarded weights over n_signatures.
        our_bytes = (n_channels * (solve_bytes + 2 * score_bytes)
                     + n_signatures * (1 + solve_bytes) + 16)
        # cuML's internal per-problem workspace (estimate): active submatrix
        # (n_channels x n_sig) plus normal-equation scratch (~n_sig^2), in the solve dtype.
        cuml_bytes = (n_channels * n_signatures + n_signatures * n_signatures) * solve_bytes
        per_problem = our_bytes + _CHUNK_WORKSPACE_FACTOR * cuml_bytes
        budget = _CHUNK_MEM_FRACTION * _GPU_POOL_BYTES
        chunk = int(budget // per_problem) if per_problem > 0 else _CHUNK_MAX
        chunk = int(min(_CHUNK_MAX, max(_CHUNK_MIN, chunk)))
        self.auto_batch_size = chunk
        if not getattr(self, '_chunk_logged', False):  # log once per process
            print(f"GPU auto batch size: {chunk} problems/call "
                  f"(n_channels={n_channels}, n_signatures={n_signatures}, "
                  f"pool={_GPU_POOL_BYTES / (1 << 30):.1f} GiB)")
            self._chunk_logged = True
        return chunk

    def solve_and_score(self, A, B, norm_obs, masks, b_index, chunk_size):
        import cupy as cp
        A_dev, B_solve, B_score, norm_dev = self._resident(A, B, norm_obs)

        n_problems = masks.shape[1]
        sims = np.empty(n_problems)
        for start in range(0, n_problems, chunk_size):
            stop = min(start + chunk_size, n_problems)
            # Our masks are C-contiguous, so force F-contiguity on the device (small,
            # boolean) - cuML wants column-major masks and reads a block per problem.
            masks_dev = cp.asfortranarray(cp.asarray(masks[:, start:stop], dtype=cp.uint8))
            bi_dev = cp.asarray(b_index[start:stop])
            _, fitted = self._nnls_batched(A_dev, B_solve, masks_dev,
                                           b_index=bi_dev, compute_fitted=True)
            if self._mixed:
                fitted = fitted.astype(self._score_dtype, copy=False)
            resid = cp.linalg.norm(B_score[:, bi_dev] - fitted, axis=0)
            sims[start:stop] = cp.asnumpy(1.0 - resid / norm_dev[bi_dev])
        return sims


def configure_gpu_backend(precision='mixed'):
    """Enable the cuML GPU backend for batched_solve_and_score().

    Called once from main() when --use_gpu is set. Construction is deferred to here
    (not import time) so that importing this module, and the whole CPU path, never
    require cupy/cuml to be installed. See _GpuBatchedSolver for the precision modes.
    """
    global _GPU_SOLVER
    _GPU_SOLVER = _GpuBatchedSolver(precision=precision)


def batched_solve_and_score(A, B, norm_obs, masks, b_index, chunk_size=4096):
    """Solve a batch of masked NNLS problems and return each fit's similarity.

    THIS IS THE GPU SWAP POINT. It is the whole of the optimisation hot path:
    ~99.9% of all NNLS solves in a run happen inside this call, and it is the only
    function a GPU/cuML backend needs to implement.

    For each problem j it solves NNLS(A[:, masks[:, j]], B[:, b_index[j]]) and
    returns a single scalar, the L2_normalised_by_first similarity of the fit:

        sims[j] = 1 - ||B[:, b_index[j]] - fitted_j|| / norm_obs[b_index[j]]

    The weights and the fitted vector are *internal* - the greedy caller never
    sees them. A backend should therefore fuse the solve with the residual-norm
    reduction and keep both on the device, returning only ``sims``. Handing back
    ``fitted`` instead would move n_channels * 8 bytes per problem to deliver 8
    bytes of signal: an amplification of exactly n_channels (96x at SBS-96, 4608x
    at SBS-4608), independent of the number of signatures. It would also force the
    host to reduce all of that data. Both costs grow with n_signatures, since the
    greedy issues ~k*(k+1)/2 problems per sample.

    ``A``, ``B`` and ``norm_obs`` are loop-invariant across the greedy steps of a
    given batch, so a backend is free to keep them device-resident between calls;
    only ``masks``/``b_index`` change per call, and only ``sims`` comes back.

    Parameters
    ----------
    A : ndarray (n_channels, n_signatures)
        Full signature matrix, shared by all problems.
    B : ndarray (n_channels, n_targets)
        Distinct target spectra; not duplicated per problem.
    norm_obs : ndarray (n_targets,)
        ||b|| for each target, precomputed.
    masks : ndarray (n_signatures, n_problems), bool
        Column j selects the active signatures of problem j.
    b_index : ndarray (n_problems,), int
        Target column of B for each problem.
    chunk_size : int
        CPU implementation detail: bounds the intermediate (n_channels x chunk)
        buffers. A GPU backend may ignore it and choose its own tiling.

    Returns
    -------
    ndarray (n_problems,) : similarity of each fit to its target.
    """
    if _GPU_SOLVER is not None:
        # cuML backend: solve + fused residual-norm reduction on the device, with
        # A/B/norm_obs kept resident across greedy steps and the batch tiled by
        # chunk_size (device-memory bound, as on the CPU path).
        return _GPU_SOLVER.solve_and_score(A, B, norm_obs, masks, b_index, chunk_size)

    n_problems = masks.shape[1]
    sims = np.empty(n_problems)
    for start in range(0, n_problems, chunk_size):
        stop = min(start + chunk_size, n_problems)
        bi = b_index[start:stop]
        _, fitted = nnls_batched(A, B, masks[:, start:stop], b_index=bi)
        sims[start:stop] = _l2_similarity_batched(B[:, bi], fitted, norm_obs[bi])
    return sims


def remove_weak_signatures_batched(sig_values, B, norm_obs, masks,
                                   weak_threshold=0.01, idx=8, chunk_size=4096):
    """Sample-transposed greedy leave-one-out removal for a whole batch of samples.

    Mathematically identical to running remove_weak_signatures() on each sample
    independently: every sample keeps its own mask, its own base similarity and
    its own stopping decision, and no quantity is ever pooled across samples.
    Samples drop out of the active set as they converge, so the batch shrinks.

    Each greedy step issues one batched_solve_and_score() call covering the
    leave-one-out trials of *all* still-active samples - this is where essentially
    all of the run's NNLS work happens, and the only place a GPU backend is needed.

    Parameters
    ----------
    sig_values : (n_channels, n_signatures)
    B : (n_channels, N) observed spectra, one column per sample
    norm_obs : (N,) precomputed ||b|| per sample
    masks : (n_signatures, N) bool, initial active set per sample
    chunk_size : int, max NNLS problems solved per batched call (memory bound)

    Returns
    -------
    masks : (n_signatures, N) bool, final active set per sample
    """
    if idx != 8:
        raise ValueError("Batched removal only supports the L2_normalised_by_first metric")

    masks = masks.copy()
    n_sig, N = masks.shape

    # Initial base similarity for every sample (fit on its starting mask)
    base_sim = batched_solve_and_score(sig_values, B, norm_obs, masks, np.arange(N), chunk_size)

    # A sample stops once it is down to a single signature
    active = np.flatnonzero(masks.sum(axis=0) > 1)

    # Group samples so that one trial batch stays within chunk_size problems
    samples_per_group = max(1, chunk_size // max(1, n_sig))

    while active.size:
        still_active = []
        for group_start in range(0, active.size, samples_per_group):
            group = active[group_start:group_start + samples_per_group]

            cols_per_sample = [np.flatnonzero(masks[:, j]) for j in group]
            counts = np.array([c.size for c in cols_per_sample])
            offsets = np.zeros(group.size + 1, dtype=np.int64)
            np.cumsum(counts, out=offsets[1:])
            n_trials = int(offsets[-1])

            # Build the leave-one-out trials: trial p of sample j drops its p-th active column
            trial_masks = np.empty((n_sig, n_trials), dtype=bool)
            trial_bidx = np.empty(n_trials, dtype=np.int64)
            for s, j in enumerate(group):
                lo, hi = offsets[s], offsets[s + 1]
                trial_masks[:, lo:hi] = masks[:, j][:, None]
                trial_masks[cols_per_sample[s], np.arange(lo, hi)] = False
                trial_bidx[lo:hi] = j

            sims = batched_solve_and_score(sig_values, B, norm_obs, trial_masks, trial_bidx, chunk_size)

            for s, j in enumerate(group):
                lo, hi = offsets[s], offsets[s + 1]
                # Contribution of each signature = drop in similarity if removed
                contributions = base_sim[j] - sims[lo:hi]
                p_weakest = int(np.argmin(contributions))
                if contributions[p_weakest] < weak_threshold:
                    masks[cols_per_sample[s][p_weakest], j] = False
                    # The chosen trial *is* the fit on the new mask, so its
                    # similarity is exactly the recomputed base similarity.
                    base_sim[j] = sims[lo + p_weakest]
                    if masks[:, j].sum() > 1:
                        still_active.append(j)
                # else: no signature is weak enough -> this sample has converged

        active = np.asarray(still_active, dtype=np.int64)

    return masks


def optimise_signatures_batched(B, norm_obs, sig_values, initial_cols, strategy='removal',
                                weak_threshold=0.01, similarity_index=-3, chunk_size=4096):
    """Batched (sample-transposed) signature optimisation.

    Returns a list of per-sample column-index lists (ascending, matching the
    scalar path's ordering), or None if this strategy/metric combination is not
    supported by the batched path and the caller should fall back per sample.
    """
    idx = similarity_index % 11
    if strategy != 'removal' or idx != 8:
        return None

    masks = np.zeros((sig_values.shape[1], B.shape[1]), dtype=bool)
    masks[initial_cols, :] = True
    masks = remove_weak_signatures_batched(sig_values, B, norm_obs, masks,
                                           weak_threshold, idx, chunk_size)
    return [np.flatnonzero(masks[:, j]).tolist() for j in range(B.shape[1])]


def _similarity_from_fitted(observed, fitted, idx, norm_obs):
    """Return ONLY the single similarity metric (position ``idx`` in the
    stat_info layout) from an already-computed fitted vector.

    Shared by the single-solve and batched paths so they use identical
    arithmetic. The default optimisation metric, L2_normalised_by_first (idx 8),
    is computed directly from the residual without scipy.spatial.distance.
    """
    if idx == 8:  # L2_normalised_by_first: 1 - ||observed - fitted|| / ||observed||
        return 1.0 - np.linalg.norm(observed - fitted) / norm_obs
    if idx in _IDX_TO_METRIC:
        metric, normalise = _IDX_TO_METRIC[idx]
        return calculate_similarity(observed, fitted, metric=metric, normalise=normalise)
    # Fallback for any other (non-similarity) index: reproduce the exact value the
    # original stat_info list would have held at that position.
    rss, chi2, r2 = calculate_stats_numba(observed, fitted)
    stat_info = [observed.sum(), rss, chi2, r2,
                 calculate_similarity(observed, fitted),
                 calculate_similarity(observed, fitted, metric='Correlation'),
                 calculate_similarity(observed, fitted, metric='Chebyshev', normalise=True),
                 calculate_similarity(observed, fitted, metric='L1', normalise=True),
                 calculate_similarity(observed, fitted, metric='L2_normalised_by_first'),
                 calculate_similarity(observed, fitted, metric='L3', normalise=True),
                 calculate_similarity(observed, fitted, metric='jensen-shannon')]
    return stat_info[idx]


def _attribution_similarity(sig_sub, observed, idx, norm_obs):
    """Solve a single NNLS problem and return only the optimisation metric.

    Used by the scalar (per-sample) optimisation path for its base/final solves.
    The GPU path scores whole batches through batched_solve_and_score() instead.
    """
    weights, _ = nnls(sig_sub, observed)
    fitted = sig_sub @ weights
    return _similarity_from_fitted(observed, fitted, idx, norm_obs)


def perform_signature_attribution(selected_mutations, signatures, normalise_mutations=False,
                                  verbose=False, weights=None):
    """Optimized NNLS attribution with cached computations.

    Computes the full stat_info list (burden + all similarities). This is used
    for the final per-sample attribution that is written to disk; the iterative
    optimisation loops use _attribution_similarity() instead, which computes only
    the single metric they need.

    If ``weights`` is supplied (e.g. from a batched NNLS solve), the NNLS step is
    skipped and those weights are used directly. This lets the non-optimised path
    solve all samples through batched_nnls_shared_A() and still reuse this
    function for the per-sample post-processing.
    """
    if signatures.empty:
        if verbose:
            print('Zero signatures provided to NNLS.')
        n_stats = 11
        return np.nan, np.nan, np.nan, np.nan, [sum(selected_mutations)] + [np.nan] * (n_stats - 1)

    # Convert to numpy for speed
    sig_array = signatures.values
    mut_array = np.array(selected_mutations)

    # NNLS solve (skipped when weights are supplied by a batched solver)
    if weights is None:
        weights, _ = nnls(sig_array, mut_array)

    # Normalize weights
    weight_sum = weights.sum()
    normalised_weights = weights / weight_sum if weight_sum > 0 else weights

    mutation_numbers = (normalised_weights * mut_array.sum() if normalise_mutations
                       else weights)

    # Calculate fitted and residuals
    fitted = sig_array @ weights
    residuals = mut_array - fitted

    # Fast statistics with numba
    rss, chi2, r2 = calculate_stats_numba(mut_array, fitted)
    input_mutational_burden = mut_array.sum()

    # Similarity metrics (these are already optimized in common_methods)
    cosine_similarity = calculate_similarity(mut_array, fitted)
    correlation = calculate_similarity(mut_array, fitted, metric='Correlation')
    chebyshev_similarity = calculate_similarity(mut_array, fitted, metric='Chebyshev', normalise=True)
    L1_similarity = calculate_similarity(mut_array, fitted, metric='L1', normalise=True)
    L2_similarity = calculate_similarity(mut_array, fitted, metric='L2_normalised_by_first')
    L3_similarity = calculate_similarity(mut_array, fitted, metric='L3', normalise=True)
    jensenshannon_similarity = calculate_similarity(mut_array, fitted, metric='jensen-shannon')

    stat_info = [input_mutational_burden, rss, chi2, r2, cosine_similarity,
                 correlation, chebyshev_similarity, L1_similarity, L2_similarity,
                 L3_similarity, jensenshannon_similarity]

    if verbose:
        if chi2 > 1e10:
            print('************* High chi2 sample *************')
        print(f'Signatures: {signatures.columns.tolist()}')
        print(f'Stats: {stat_info}')

    return normalised_weights, mutation_numbers, fitted, residuals, stat_info


def remove_weak_signatures(observed, norm_obs, sig_values, active_cols, col_names,
                           weak_threshold=0.01, idx=8, verbose=False):
    """Greedy leave-one-out removal of weak signatures.

    Operates on integer column indices into ``sig_values`` (no pandas inside the
    loop) and evaluates a single similarity metric per candidate. Column order in
    ``active_cols`` is preserved, so the surviving set matches the original
    DataFrame-based implementation.

    Returns (active_cols, final_similarity).
    """
    active_cols = list(active_cols)
    base_similarity = _attribution_similarity(sig_values[:, active_cols], observed, idx, norm_obs)

    if verbose:
        print(f'Starting removal loop. Base similarity: {base_similarity}')
        print(f'Current signatures: {[col_names[c] for c in active_cols]}')

    while len(active_cols) > 1:
        # Leave-one-out candidates: each drops one signature, all sharing the same
        # target b. (In GPU mode this loop is replaced wholesale by the
        # sample-transposed remove_weak_signatures_batched(), which batches across
        # samples instead of across one sample's candidates.)
        column_sets = [active_cols[:p] + active_cols[p + 1:] for p in range(len(active_cols))]
        candidate_weights = [nnls(sig_values[:, cols], observed)[0] for cols in column_sets]

        # Contribution of each signature = drop in similarity if it is removed.
        contributions = np.empty(len(active_cols))
        for p, cols in enumerate(column_sets):
            fitted = sig_values[:, cols] @ candidate_weights[p]
            contributions[p] = base_similarity - _similarity_from_fitted(observed, fitted, idx, norm_obs)

        p_weakest = int(np.argmin(contributions))
        if contributions[p_weakest] < weak_threshold:
            dropped = active_cols.pop(p_weakest)
            base_similarity = _attribution_similarity(sig_values[:, active_cols], observed, idx, norm_obs)
            if verbose:
                print(f'Dropped {col_names[dropped]}, {len(active_cols)} signatures left')
        else:
            if verbose:
                print('All weak signatures removed')
            break

    if len(active_cols) <= 1 and verbose:
        print(f'Only one signature left: {[col_names[c] for c in active_cols]}')

    final_similarity = _attribution_similarity(sig_values[:, active_cols], observed, idx, norm_obs)
    return active_cols, final_similarity


def remove_weak_signatures_trajectory(observed, norm_obs, sig_values, active_cols, idx=8):
    """Greedy leave-one-out removal recorded to full depth, penalty-free.

    Runs exactly the same removal loop as remove_weak_signatures() but with no
    stopping threshold: it removes the weakest signature at every step, down to a
    single one, recording the contribution (similarity drop) at which each is
    removed. Because the removal *choice* at each step (argmin contribution) does not
    depend on the penalty - only the stop/continue decision does - a removal run at
    any weak_threshold is a prefix of this one trajectory. select_from_trajectory()
    can therefore reconstruct the surviving set at ANY penalty from this single pass,
    exactly reproducing remove_weak_signatures() without re-solving per penalty.

    Returns
    -------
    active_cols : list[int]
        the starting set, in original order.
    trajectory : list[(int, float)]
        (dropped_col, contribution) in removal order, length len(active_cols) - 1.
    """
    active_cols = list(active_cols)
    remaining = list(active_cols)
    trajectory = []
    base_similarity = _attribution_similarity(sig_values[:, remaining], observed, idx, norm_obs)

    while len(remaining) > 1:
        column_sets = [remaining[:p] + remaining[p + 1:] for p in range(len(remaining))]
        contributions = np.empty(len(remaining))
        for p, cols in enumerate(column_sets):
            fitted = sig_values[:, cols] @ nnls(sig_values[:, cols], observed)[0]
            contributions[p] = base_similarity - _similarity_from_fitted(observed, fitted, idx, norm_obs)
        p_weakest = int(np.argmin(contributions))
        trajectory.append((remaining.pop(p_weakest), float(contributions[p_weakest])))
        base_similarity = _attribution_similarity(sig_values[:, remaining], observed, idx, norm_obs)

    return active_cols, trajectory


def select_from_trajectory(active_cols, trajectory, weak_threshold):
    """Reconstruct the surviving signature set at ``weak_threshold`` from a trajectory.

    Replays the recorded removals: drop each signature while its recorded
    contribution is < weak_threshold, and stop at the first contribution that is not
    (mirroring remove_weak_signatures()'s break). Returns the survivors in the
    original ``active_cols`` order. This is a cheap O(k) re-threshold, so the optimal
    penalty can be scanned at arbitrary resolution from one trajectory per sample.
    """
    removed = set()
    for dropped, contribution in trajectory:
        if contribution < weak_threshold:
            removed.add(dropped)
        else:
            break
    return [c for c in active_cols if c not in removed]


def remove_weak_signatures_batched_trajectory(sig_values, B, norm_obs, masks, idx=8, chunk_size=4096):
    """Batched (GPU-path) counterpart of remove_weak_signatures_trajectory().

    Runs exactly the greedy steps of remove_weak_signatures_batched() - same trial
    batches, same base-similarity update from the chosen trial - but never stops: every
    sample is followed down to a single signature. select_from_trajectory() on the
    result reproduces remove_weak_signatures_batched() at any weak_threshold, provided
    the solver's per-problem results do not depend on how problems are batched
    (checked on real GPUs by tests/gpu_acceptance.sh).

    Returns
    -------
    list of trajectories, one per column of B: [(dropped_col, contribution), ...]
    """
    if idx != 8:
        raise ValueError("Batched removal only supports the L2_normalised_by_first metric")

    masks = masks.copy()
    n_sig, N = masks.shape
    base_sim = batched_solve_and_score(sig_values, B, norm_obs, masks, np.arange(N), chunk_size)
    trajectories = [[] for _ in range(N)]
    active = np.flatnonzero(masks.sum(axis=0) > 1)
    samples_per_group = max(1, chunk_size // max(1, n_sig))

    while active.size:
        still_active = []
        for group_start in range(0, active.size, samples_per_group):
            group = active[group_start:group_start + samples_per_group]

            cols_per_sample = [np.flatnonzero(masks[:, j]) for j in group]
            counts = np.array([c.size for c in cols_per_sample])
            offsets = np.zeros(group.size + 1, dtype=np.int64)
            np.cumsum(counts, out=offsets[1:])
            n_trials = int(offsets[-1])

            trial_masks = np.empty((n_sig, n_trials), dtype=bool)
            trial_bidx = np.empty(n_trials, dtype=np.int64)
            for s, j in enumerate(group):
                lo, hi = offsets[s], offsets[s + 1]
                trial_masks[:, lo:hi] = masks[:, j][:, None]
                trial_masks[cols_per_sample[s], np.arange(lo, hi)] = False
                trial_bidx[lo:hi] = j

            sims = batched_solve_and_score(sig_values, B, norm_obs, trial_masks, trial_bidx, chunk_size)

            for s, j in enumerate(group):
                lo, hi = offsets[s], offsets[s + 1]
                contributions = base_sim[j] - sims[lo:hi]
                p_weakest = int(np.argmin(contributions))
                dropped = int(cols_per_sample[s][p_weakest])
                trajectories[j].append((dropped, float(contributions[p_weakest])))
                masks[dropped, j] = False
                base_sim[j] = sims[lo + p_weakest]
                if masks[:, j].sum() > 1:
                    still_active.append(j)

        active = np.asarray(still_active, dtype=np.int64)

    return trajectories


# Percentiles of the pooled trajectory breakpoints spanned by the automatic thresholds
AUTO_THRESHOLD_PERCENTILES = (5, 95)
# Breakpoints below this are numerical noise (removing signatures with zero weight, whose
# contribution is +-1e-16 in fp64 and ~1e-8 on a GPU in fp32) and are ignored; the
# no-removal end of the scan is covered by always including 0.
AUTO_THRESHOLD_FLOOR = 1e-6


def auto_weak_thresholds(trajectories, n_thresholds=12, percentiles=AUTO_THRESHOLD_PERCENTILES,
                         floor=AUTO_THRESHOLD_FLOOR):
    """Weak thresholds to scan, chosen from removal trajectories.

    A sample's surviving signature set only changes where the threshold crosses a
    running maximum of its removal contributions (select_from_trajectory() stops at the
    first contribution >= threshold). These breakpoints are pooled over all samples;
    the result is 0 (plain NNLS) plus n_thresholds values log-spaced between the given
    percentiles of the breakpoints >= floor, rounded to 2 significant digits.
    """
    breakpoints = [np.maximum.accumulate([c for _, c in t]) for t in trajectories if t]
    breakpoints = np.concatenate(breakpoints) if breakpoints else np.empty(0)
    breakpoints = breakpoints[breakpoints >= floor]
    if breakpoints.size == 0 or n_thresholds < 1:
        return [0.0]
    low, high = np.percentile(breakpoints, percentiles)
    values = [low] if n_thresholds == 1 or high <= low else np.geomspace(low, high, n_thresholds)
    return [0.0] + sorted({float(f'{v:.2g}') for v in values})


def removal_trajectories(input_mutations, signatures, sel_sig_nums, args):
    """One removal trajectory per sample, {sample position: (active_cols, trajectory)}.

    Mirrors process_samples_batch(): the batched trajectory in GPU mode, the per-sample
    one otherwise, with the same inputs, so that process_samples_batch(...,
    trajectories=...) reproduces a per-threshold run at every weak threshold.
    Samples without mutations are left out, as process_samples_batch() skips them.
    """
    initial_signatures = signatures.iloc[:, sel_sig_nums]
    name_to_col = {name: i for i, name in enumerate(signatures.columns)}
    sig_values = np.asarray(signatures.values, dtype=float)
    initial_cols = [name_to_col[name] for name in initial_signatures.columns]

    valid = [(s_i, input_mutations[sample].values) for s_i, sample in enumerate(input_mutations.columns)
             if not input_mutations[sample].values.sum() <= 0]  # same test as process_samples_batch
    if not valid or initial_signatures.empty:
        return {}

    if getattr(args, 'use_gpu', False):
        B = np.column_stack([np.asarray(values, dtype=float) for _, values in valid])
        norm_obs = np.array([np.linalg.norm(B[:, j]) for j in range(B.shape[1])])
        gpu_chunk = getattr(args, 'gpu_batch_size', 0)
        if gpu_chunk <= 0:
            gpu_chunk = _GPU_SOLVER.suggest_chunk(sig_values.shape[0], sig_values.shape[1])
        masks = np.zeros((sig_values.shape[1], B.shape[1]), dtype=bool)
        masks[initial_cols, :] = True
        batched = remove_weak_signatures_batched_trajectory(sig_values, B, norm_obs, masks,
                                                            chunk_size=gpu_chunk)
        # the batched path reports surviving columns in ascending order (np.flatnonzero)
        return {s_i: (sorted(initial_cols), batched[j]) for j, (s_i, _) in enumerate(valid)}

    trajectories = {}
    for s_i, values in valid:
        observed = np.asarray(values, dtype=float)
        trajectories[s_i] = remove_weak_signatures_trajectory(observed, np.linalg.norm(observed),
                                                              sig_values, initial_cols)
    return trajectories


def add_strong_signatures(observed, norm_obs, sig_values, active_cols, candidate_cols, col_names,
                          strong_threshold=0.05, idx=8, verbose=False):
    """Greedy addition of strong signatures drawn from ``candidate_cols``.

    Like remove_weak_signatures, this works on integer column indices and a
    single metric. The surviving set is returned sorted by signature name, to
    match the original implementation's reindex(sorted(columns)).

    Returns (active_cols, final_similarity).
    """
    active_cols = list(active_cols)
    active_set = set(active_cols)
    # Candidates not already in the model, keeping the original column ordering.
    remaining = [c for c in candidate_cols if c not in active_set]

    if not remaining:
        sim = _attribution_similarity(sig_values[:, active_cols], observed, idx, norm_obs)
        return active_cols, sim

    base_similarity = _attribution_similarity(sig_values[:, active_cols], observed, idx, norm_obs)
    if verbose:
        print(f'Starting addition loop. Base similarity: {base_similarity}')

    while remaining:
        # Add-one candidates: each adds one signature, all sharing the same target b.
        column_sets = [active_cols + [c] for c in remaining]
        candidate_weights = [nnls(sig_values[:, cols], observed)[0] for cols in column_sets]

        contributions = np.empty(len(remaining))
        for p, cols in enumerate(column_sets):
            fitted = sig_values[:, cols] @ candidate_weights[p]
            contributions[p] = _similarity_from_fitted(observed, fitted, idx, norm_obs) - base_similarity

        p_strongest = int(np.argmax(contributions))
        if contributions[p_strongest] > strong_threshold:
            chosen = remaining.pop(p_strongest)
            active_cols.append(chosen)
            base_similarity = _attribution_similarity(sig_values[:, active_cols], observed, idx, norm_obs)
            if verbose:
                print(f'Added {col_names[chosen]}')
        else:
            if verbose:
                print('All strong signatures added')
            break

    # Match the original: reindex(sorted(columns)) -> order surviving sigs by name.
    active_cols = sorted(active_cols, key=lambda c: col_names[c])
    final_similarity = _attribution_similarity(sig_values[:, active_cols], observed, idx, norm_obs)
    return active_cols, final_similarity


def optimise_signatures(selected_mutations, initial_signatures, all_available_signatures,
                       strategy='removal', weak_threshold=0.01, strong_threshold=0.05,
                       similarity_index=-3, loops_limit=100, verbose=False):
    """Optimise the signature set for one sample.

    The candidate search runs entirely on numpy arrays / integer column indices
    and evaluates a single similarity metric per NNLS solve (see
    _attribution_similarity). The result is rebuilt as a DataFrame once, with the
    same column ordering as the original implementation, so downstream output is
    unchanged.
    """
    observed = np.asarray(selected_mutations, dtype=float)
    norm_obs = np.linalg.norm(observed)
    idx = similarity_index % 11  # default -3 maps to 8 (L2_normalised_by_first)

    all_names = list(all_available_signatures.columns)
    name_to_col = {name: i for i, name in enumerate(all_names)}
    sig_values = np.asarray(all_available_signatures.values, dtype=float)
    candidate_cols = list(range(len(all_names)))

    active_cols = [name_to_col[name] for name in initial_signatures.columns]

    if not active_cols:
        return all_available_signatures.iloc[:, :0]

    if verbose:
        base_similarity = _attribution_similarity(sig_values[:, active_cols], observed, idx, norm_obs)
        print(f'Initial signatures: {[all_names[c] for c in active_cols]}')
        print(f'Initial similarity: {base_similarity}')

    if strategy == 'removal':
        active_cols, _ = remove_weak_signatures(
            observed, norm_obs, sig_values, active_cols, all_names,
            weak_threshold, idx, verbose)

    elif strategy == 'addition':
        active_cols, _ = add_strong_signatures(
            observed, norm_obs, sig_values, active_cols, candidate_cols, all_names,
            strong_threshold, idx, verbose)

    elif strategy == 'add-remove':
        converging_similarity = _attribution_similarity(
            sig_values[:, active_cols], observed, idx, norm_obs)

        for loop_counter in range(1, loops_limit + 1):
            prev_similarity = converging_similarity

            active_cols, converging_similarity = add_strong_signatures(
                observed, norm_obs, sig_values, active_cols, candidate_cols, all_names,
                strong_threshold, idx, verbose)

            active_cols, converging_similarity = remove_weak_signatures(
                observed, norm_obs, sig_values, active_cols, all_names,
                weak_threshold, idx, verbose)

            convergence_delta = abs(converging_similarity - prev_similarity)

            if verbose:
                print(f'Loop {loop_counter} done. Convergence delta: {convergence_delta}')

            if convergence_delta <= 0:
                break

            if loop_counter >= loops_limit:
                warnings.warn(f"Maximum iterations ({loops_limit}) reached")
                break
    else:
        raise ValueError(f'Unknown strategy: {strategy}')

    final_names = [all_names[c] for c in active_cols]
    if verbose:
        print(f'Final signatures: {final_names}')

    return all_available_signatures.loc[:, final_names]


def process_samples_batch(input_mutations, signatures, sel_sig_nums, args, trajectories=None):
    """Process all samples, accumulating results in preallocated numpy arrays.

    Per-sample results are written into numpy arrays and assembled into output
    DataFrames once at the end. The original implementation assigned into
    DataFrames row-by-row via .loc, which is slow at scale.

    ``trajectories`` (from removal_trajectories(), removal strategy only) replaces the
    greedy optimisation by select_from_trajectory() at args.weak_threshold; everything
    else is unchanged, so the output is identical to a run without it.
    """
    samples = list(input_mutations.columns)
    signature_columns_list = signatures.columns[sel_sig_nums].tolist()
    col_pos = {name: i for i, name in enumerate(signature_columns_list)}

    n_samples = len(samples)
    n_sigs = len(signature_columns_list)
    n_channels = len(input_mutations.index)

    # Pre-allocate output arrays (rows kept zero for skipped/empty samples)
    weights_arr = np.zeros((n_samples, n_sigs))
    mutations_arr = np.zeros((n_samples, n_sigs))
    stat_arr = np.zeros((n_samples, 11))
    residuals_arr = np.zeros((n_channels, n_samples))
    fitted_arr = np.zeros((n_channels, n_samples))

    initial_signatures = signatures.iloc[:, sel_sig_nums]

    def _store(s_i, cols, normalised_weights, mutation_numbers, stat_info, residuals, fitted):
        weights_arr[s_i, cols] = normalised_weights
        mutations_arr[s_i, cols] = mutation_numbers
        stat_arr[s_i, :] = stat_info
        residuals_arr[:, s_i] = residuals
        fitted_arr[:, s_i] = fitted

    # GPU mode: run the greedy optimisation transposed across samples, so each
    # batched NNLS call covers the leave-one-out trials of every active sample.
    use_gpu_batched = (args.optimise_signatures
                       and getattr(args, 'use_gpu', False)
                       and args.optimisation_strategy == 'removal'
                       and not initial_signatures.empty)

    if not args.optimise_signatures and not initial_signatures.empty:
        # Non-optimised: every sample shares the same signature matrix, so all
        # NNLS solves go through a single batched call (the shared-A GPU swap
        # point). Per-sample post-processing reuses perform_signature_attribution
        # via its weights= argument, keeping output identical to the per-sample path.
        sig_array = np.asarray(initial_signatures.values, dtype=float)
        cols = [col_pos[name] for name in initial_signatures.columns]

        valid_s_i, valid_b = [], []
        for s_i, sample in enumerate(samples):
            selected_mutations = input_mutations[sample].values
            if selected_mutations.sum() <= 0:
                warnings.warn(f"Sample {sample}: Zero mutations, skipping")
                continue
            valid_s_i.append(s_i)
            valid_b.append(np.asarray(selected_mutations, dtype=float))

        if valid_b:
            B = np.column_stack(valid_b)
            weights_batch = batched_nnls_shared_A(sig_array, B)
            for j, s_i in enumerate(valid_s_i):
                normalised_weights, mutation_numbers, fitted, residuals, stat_info = \
                    perform_signature_attribution(B[:, j], initial_signatures,
                                                normalise_mutations=args.normalise_mutations,
                                                verbose=args.verbose, weights=weights_batch[:, j])
                _store(s_i, cols, normalised_weights, mutation_numbers, stat_info, residuals, fitted)

    elif use_gpu_batched:
        all_names = list(signatures.columns)
        name_to_col = {name: i for i, name in enumerate(all_names)}
        sig_values = np.asarray(signatures.values, dtype=float)
        initial_cols = [name_to_col[name] for name in initial_signatures.columns]

        valid_s_i, valid_b = [], []
        for s_i, sample in enumerate(samples):
            selected_mutations = input_mutations[sample].values
            if selected_mutations.sum() <= 0:
                warnings.warn(f"Sample {sample}: Zero mutations, skipping")
                continue
            valid_s_i.append(s_i)
            valid_b.append(np.asarray(selected_mutations, dtype=float))

        if valid_b:
            B = np.column_stack(valid_b)
            # match the scalar path's per-sample norm computation exactly
            norm_obs = np.array([np.linalg.norm(B[:, j]) for j in range(B.shape[1])])

            if trajectories is not None:
                per_sample_cols = [select_from_trajectory(*trajectories[s_i], args.weak_threshold)
                                   for s_i in valid_s_i]
            else:
                gpu_chunk = getattr(args, 'gpu_batch_size', 0)
                if gpu_chunk <= 0:  # auto-size the batch to the GPU's memory
                    gpu_chunk = _GPU_SOLVER.suggest_chunk(sig_values.shape[0], sig_values.shape[1])
                per_sample_cols = optimise_signatures_batched(
                    B, norm_obs, sig_values, initial_cols,
                    strategy=args.optimisation_strategy,
                    weak_threshold=args.weak_threshold,
                    chunk_size=gpu_chunk)

            for j, s_i in enumerate(valid_s_i):
                final_names = [all_names[c] for c in per_sample_cols[j]]
                if not final_names:
                    continue
                final_signatures = signatures.loc[:, final_names]

                normalised_weights, mutation_numbers, fitted, residuals, stat_info = \
                    perform_signature_attribution(B[:, j], final_signatures,
                                                normalise_mutations=args.normalise_mutations,
                                                verbose=args.verbose)

                cols = [col_pos[name] for name in final_names]
                _store(s_i, cols, normalised_weights, mutation_numbers, stat_info, residuals, fitted)
    else:
        for s_i, sample in enumerate(samples):
            selected_mutations = input_mutations[sample].values

            if selected_mutations.sum() <= 0:
                warnings.warn(f"Sample {sample}: Zero mutations, skipping")
                continue

            if args.optimise_signatures and trajectories is not None:
                keep = select_from_trajectory(*trajectories[s_i], args.weak_threshold)
                final_signatures = signatures.loc[:, [signatures.columns[c] for c in keep]]
            elif args.optimise_signatures:
                final_signatures = optimise_signatures(
                    selected_mutations, initial_signatures, signatures,
                    strategy=args.optimisation_strategy,
                    weak_threshold=args.weak_threshold,
                    strong_threshold=args.strong_threshold,
                    verbose=args.verbose)
            else:
                final_signatures = initial_signatures

            if not final_signatures.empty:
                normalised_weights, mutation_numbers, fitted, residuals, stat_info = \
                    perform_signature_attribution(selected_mutations, final_signatures,
                                                normalise_mutations=args.normalise_mutations,
                                                verbose=args.verbose)

                cols = [col_pos[name] for name in final_signatures.columns]
                _store(s_i, cols, normalised_weights, mutation_numbers, stat_info, residuals, fitted)

    stat_columns = ['Mutational burden', 'RSS', 'Chi2', 'R2',
                    'Cosine similarity', 'Correlation', 'Chebyshev similarity',
                    'L1 similarity', 'L2 similarity', 'L3 similarity',
                    'Jensen-Shannon similarity']

    output_weights = pd.DataFrame(weights_arr, index=samples, columns=signature_columns_list)
    output_mutations = pd.DataFrame(mutations_arr, index=samples, columns=signature_columns_list)
    output_stat_info = pd.DataFrame(stat_arr, index=samples, columns=stat_columns)
    residuals_dataframe = pd.DataFrame(residuals_arr, index=input_mutations.index, columns=samples)
    fitted_dataframe = pd.DataFrame(fitted_arr, index=input_mutations.index, columns=samples)

    output_weights.index.name = output_mutations.index.name = output_stat_info.index.name = 'Sample'

    return output_weights, output_mutations, output_stat_info, residuals_dataframe, fitted_dataframe


def attribute_at_thresholds(input_mutations, signatures, sel_sig_nums, args,
                            weak_thresholds, strong_thresholds, trajectories=None):
    """process_samples_batch() outputs for every (weak, strong) threshold pair.

    With the removal strategy and several weak thresholds, the greedy removal runs once
    per sample (removal_trajectories(), unless args.threshold_scan == 'separate') and each
    threshold replays it; other strategies run once per pair as before. The removal
    strategy ignores the strong threshold, so its results are shared across strong values.
    Returns ({(weak, strong): outputs}, trajectories or None).
    """
    removal = args.optimise_signatures and args.optimisation_strategy == 'removal'
    if (trajectories is None and removal and len(weak_thresholds) > 1
            and getattr(args, 'threshold_scan', 'trajectory') == 'trajectory'):
        trajectories = removal_trajectories(input_mutations, signatures, sel_sig_nums, args)

    results = {}
    for weak in weak_thresholds:
        for strong in strong_thresholds:
            if removal and (weak, strong_thresholds[0]) in results:
                results[(weak, strong)] = results[(weak, strong_thresholds[0])]
                continue
            pair_args = argparse.Namespace(**vars(args))
            pair_args.weak_threshold, pair_args.strong_threshold = weak, strong
            results[(weak, strong)] = process_samples_batch(input_mutations, signatures, sel_sig_nums,
                                                            pair_args, trajectories)
    return results, trajectories


def bootstrap_attributions(input_mutations, signatures, sel_sig_nums, args,
                           weak_thresholds, strong_thresholds, central=None):
    """attribute_at_thresholds() outputs of one bootstrap resample of input_mutations
    (args.bootstrap_method), by (weak, strong) threshold pair. Seed the global RNG first.

    For bootstrap_residuals, central holds attribute_at_thresholds() of the original table
    (aligned to the signatures): each pair resamples the residuals of its own central fit,
    with one channel resample shared by all pairs.
    """
    method = args.bootstrap_method
    if method != 'bootstrap_residuals':
        resampled = bootstrap_mutation_table(input_mutations, method=method).reindex(signatures.index)
        return attribute_at_thresholds(resampled, signatures, sel_sig_nums, args,
                                       weak_thresholds, strong_thresholds)[0]
    if central is None:
        raise ValueError(f"{method} requires the central fit")
    rows = resample_channels(len(signatures.index))
    results, done = {}, {}
    for (weak, strong), outputs in central.items():
        # pairs sharing a central fit (removal ignores the strong threshold) share the result
        if id(outputs) not in done:
            _, _, _, residuals, fitted = outputs
            resampled = bootstrap_mutation_table(None, method, fitted=fitted, residuals=residuals, rows=rows)
            done[id(outputs)] = attribute_at_thresholds(resampled, signatures, sel_sig_nums, args,
                                                        [weak], [strong])[0][(weak, strong)]
        results[(weak, strong)] = done[id(outputs)]
    return results


def parse_thresholds(values, allow_auto=False):
    """-W/-S values -> sorted unique floats, or 'auto'."""
    if [v.lower() for v in values] == ['auto']:
        if not allow_auto:
            raise ValueError("'auto' is only supported for the weak threshold")
        return 'auto'
    thresholds = {float(v) for v in values}
    if any(np.isnan(t) for t in thresholds):
        # every comparison with NaN is false: the optimisation would silently do nothing
        raise ValueError(f"thresholds must be numbers, got {' '.join(values)}")
    return sorted(thresholds)


def write_threshold_labels(folder, weak_thresholds, strong_thresholds):
    """Write the threshold labels used in output names, one per line (for the pipeline)."""
    make_folder_if_not_exists(folder)
    for name, values in (('weak', weak_thresholds), ('strong', strong_thresholds)):
        with open(f'{folder}/{name}_thresholds.txt', 'w') as f:
            f.write(''.join(f'{format_threshold(v)}\n' for v in values))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dataset", dest="dataset_name", default='SIM')
    parser.add_argument("-t", "--mutation_type", dest="mutation_type", default='')
    parser.add_argument("-c", "--context", dest="context", default=96, type=int)
    parser.add_argument("-i", "--input_path", dest="input_path", default='input_mutation_tables/')
    parser.add_argument("-s", "--signature_path", dest="signature_tables_path", default='signature_tables/')
    parser.add_argument("-p", "--signature_prefix", dest="signatures_prefix", default='sigProfiler')
    parser.add_argument("-o", "--output_path", dest="output_path", default='output_tables/')
    parser.add_argument("-x", "--optimise_signatures", dest="optimise_signatures", action="store_true")
    parser.add_argument("-W", "--weak_threshold", dest="weak_threshold", nargs='+', default=['0.01'],
                        help="Weak threshold(s). Several values, or 'auto' (removal strategy, central run), "
                             "write one output folder per threshold (requires --add_suffix); see --threshold_scan")
    parser.add_argument("-S", "--strong_threshold", dest="strong_threshold", nargs='+', default=['0.01'],
                        help="Strong threshold(s); several values write one output folder per pair")
    parser.add_argument("-N", "--normalise_mutations", dest="normalise_mutations", action="store_true")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true")
    parser.add_argument("-n", "--number", dest="number_of_samples", default=-1, type=int)
    parser.add_argument("-B", "--bootstrap", dest="bootstrap", action="store_true")
    parser.add_argument("--bootstrap_method", dest="bootstrap_method", default='binomial',
                        choices=BOOTSTRAP_METHODS,
                        help="bootstrap_residuals resamples the residuals of the central fit, which "
                             "is recomputed in-process at each threshold pair")
    parser.add_argument("--add_suffix", dest="add_suffix", action="store_true")
    parser.add_argument("--optimisation_strategy", dest="optimisation_strategy", default='removal')
    parser.add_argument("--n_bootstrap", dest="n_bootstrap", default=0, type=int,
                       help="If > 0, run this many bootstrap iterations in a single process, "
                            "writing one indexed set of output files per iteration into "
                            "<output_path>/bootstrap_output/. Use n=1 (with --bootstrap_start_index) "
                            "for one process per iteration (CPU fan-out), or n=N for all iterations "
                            "in one process (GPU: a single reused context).")
    parser.add_argument("--bootstrap_start_index", dest="bootstrap_start_index", default=1, type=int,
                       help="Index of the first bootstrap iteration written by this process "
                            "(iterations are numbered start .. start + n_bootstrap - 1). Lets a "
                            "fanned-out / chunked launch write a distinct slice of the indices.")
    parser.add_argument("--seed", dest="seed", default=None, type=int,
                       help="Base RNG seed for bootstrap resampling. When set, iteration i is "
                            "seeded from (seed, i), so the resample for a given iteration is "
                            "identical regardless of how iterations are split across processes "
                            "(CPU fan-out vs GPU single-process). Default: non-deterministic.")
    parser.add_argument("--bootstrap_output_suffix", dest="bootstrap_output_suffix", default='',
                       help="String inserted into per-iteration bootstrap output filenames "
                            "(e.g. '<weak>_<strong>' thresholds), to match downstream expectations.")
    parser.add_argument("--use_gpu", dest="use_gpu", action="store_true",
                       help="Use the sample-transposed batched optimisation path (removal strategy), "
                            "which issues one large masked NNLS batch per greedy step across all "
                            "still-active samples. Intended for GPU offloading via "
                            "batched_solve_and_score().")
    parser.add_argument("--gpu_batch_size", dest="gpu_batch_size", default=0, type=int,
                       help="NNLS problems solved per batched GPU call. 0 (default) auto-sizes it "
                            "to the GPU's memory (larger batch for smaller contexts / bigger cards); "
                            "a positive value sets it manually. Ignored unless --use_gpu is set.")
    parser.add_argument("--gpu_precision", dest="gpu_precision", default='mixed',
                       choices=['mixed', 'fp32', 'fp64'],
                       help="Floating-point precision for the cuML GPU backend: 'mixed' "
                            "(default) solves in fp32 and scores in fp64; 'fp32'/'fp64' use "
                            "one dtype throughout. Ignored unless --use_gpu is set.")
    parser.add_argument("--threshold_scan", dest="threshold_scan", default='trajectory',
                       choices=['trajectory', 'separate'],
                       help="With several weak thresholds and the removal strategy: 'trajectory' "
                            "(default) runs the greedy removal once per sample to full depth and replays "
                            "it at each threshold (identical results); 'separate' runs it per threshold.")
    parser.add_argument("--auto_thresholds", dest="auto_thresholds", default=12, type=int,
                       help="Number of non-zero weak thresholds chosen by -W auto (plus 0).")
    parser.add_argument("--thresholds_dir", dest="thresholds_dir", default=None,
                       help="Write the threshold labels used in output names to "
                            "<dir>/weak_thresholds.txt and <dir>/strong_thresholds.txt.")
    parser.add_argument("--resolve_thresholds_only", dest="resolve_thresholds_only", action="store_true",
                       help="Only resolve the thresholds (-W auto: from the trajectories of the loaded "
                            "samples, see -n) and write them to --thresholds_dir; no attribution output.")

    args = parser.parse_args()

    try:
        weak_thresholds = parse_thresholds(args.weak_threshold, allow_auto=True)
        strong_thresholds = parse_thresholds(args.strong_threshold)
    except ValueError as e:
        parser.error(str(e))
    auto = weak_thresholds == 'auto'
    multi = args.optimise_signatures and (auto or len(weak_thresholds) > 1 or len(strong_thresholds) > 1)
    if args.resolve_thresholds_only and not args.thresholds_dir:
        parser.error("--resolve_thresholds_only requires --thresholds_dir")
    if auto and not (args.optimise_signatures and args.optimisation_strategy == 'removal'):
        parser.error("-W auto requires -x with the removal optimisation strategy")
    if auto and (args.bootstrap or args.n_bootstrap > 0):
        parser.error("-W auto is not available for bootstrap runs: every iteration must use the same "
                     "thresholds (resolve them first with --resolve_thresholds_only)")
    if multi and not args.add_suffix and not args.resolve_thresholds_only:
        parser.error("several thresholds (or -W auto) require --add_suffix (one output folder per threshold)")
    if args.resolve_thresholds_only and not auto:
        write_threshold_labels(args.thresholds_dir, weak_thresholds, strong_thresholds)
        sys.exit(0)
    if not auto:
        # single-threshold runs keep a scalar threshold in args, exactly as before
        args.weak_threshold, args.strong_threshold = weak_thresholds[0], strong_thresholds[0]

    # Setup paths
    dataset_name = args.dataset_name
    mutation_type = args.mutation_type
    context = args.context

    def output_dir(weak, strong):
        """Output folder of a (weak, strong) threshold pair."""
        path = f"{args.output_path}/{dataset_name}"
        if args.add_suffix:
            path += f'_{context}_NNLS'
            if args.optimise_signatures:
                path += f'_{format_threshold(weak)}_{format_threshold(strong)}'
            else:
                path += '_unoptimised'
        return path

    def bootstrap_suffix(weak, strong):
        if args.bootstrap_output_suffix and not multi:
            return f'_{args.bootstrap_output_suffix}'
        if args.add_suffix and args.optimise_signatures:
            # names downstream expect, as in the folder name: <weak label>_<strong label>
            # (also for a single threshold pair, e.g. -W auto resolving to 0 alone)
            return f'_{format_threshold(weak)}_{format_threshold(strong)}'
        return ''

    if not auto and not args.resolve_thresholds_only:
        for weak in weak_thresholds:
            for strong in strong_thresholds:
                make_folder_if_not_exists(output_dir(weak, strong))

    # Validate mutation type
    if not mutation_type:
        parser.error("Please specify mutation type using -t option (SBS, DBS, ID, SV, CNV)")
    elif mutation_type not in ['SBS', 'DBS', 'ID', 'SV', 'CNV']:
        raise ValueError(f"Unknown mutation type: {mutation_type}")

    # Load data based on mutation type
    sig_path = args.signature_tables_path
    inp_path = args.input_path
    sig_prefix = args.signatures_prefix

    if mutation_type == 'SBS':
        if context == 96:
            signatures = pd.read_csv(f'{sig_path}/{sig_prefix}_{mutation_type}_signatures.csv',
                                    sep=',', index_col=[0,1])
            input_mutations = pd.read_csv(f'{inp_path}/{dataset_name}/WGS_{dataset_name}.{context}.csv',
                                         sep=',', index_col=[0,1])
        elif context in [192, 288]:
            signatures = pd.read_csv(f'{sig_path}/{sig_prefix}_{mutation_type}_{context}_signatures.csv',
                                    sep=',', index_col=[0,1,2])
            input_mutations = pd.read_csv(f'{inp_path}/{dataset_name}/WGS_{dataset_name}.{context}.csv',
                                         sep=',', index_col=[0,1,2])
        elif context in [1536, 4608]:
            signatures = pd.read_csv(f'{sig_path}/{sig_prefix}_{mutation_type}_{context}_signatures.csv',
                                    sep=',', index_col=0)
            input_mutations = pd.read_csv(f'{inp_path}/{dataset_name}/WGS_{dataset_name}.{context}.csv',
                                         sep=',', index_col=0)
        else:
            raise ValueError(f"Context {context} is not supported")
    elif mutation_type == 'DBS':
        signatures = pd.read_csv(f'{sig_path}/{sig_prefix}_{mutation_type}_signatures.csv',
                                sep=',', index_col=0)
        input_mutations = pd.read_csv(f'{inp_path}/{dataset_name}/WGS_{dataset_name}.dinucs.csv',
                                     sep=',', index_col=0)
    elif mutation_type == 'ID':
        signatures = pd.read_csv(f'{sig_path}/{sig_prefix}_{mutation_type}_signatures.csv',
                                sep=',', index_col=0)
        input_mutations = pd.read_csv(f'{inp_path}/{dataset_name}/WGS_{dataset_name}.indels.csv',
                                     sep=',', index_col=0)
    else:  # SV and CNV
        signatures = pd.read_csv(f'{sig_path}/{sig_prefix}_{mutation_type}_signatures.csv',
                                sep=',', index_col=0)
        input_mutations = pd.read_csv(f'{inp_path}/{dataset_name}/WGS_{dataset_name}.{mutation_type}.csv',
                                     sep=',', index_col=0)

    print(f"Performing NNLS for {dataset_name} dataset, {mutation_type} mutation type.")
    if mutation_type == 'SBS':
        print(f"SBS context: {context}")
    if args.optimise_signatures:
        print(f"Optimised NNLS: weak thresholds = {weak_thresholds}, strong thresholds = {strong_thresholds}"
              + (f" ({args.threshold_scan} scan)" if multi else ""))
    if args.bootstrap:
        print(f"Bootstrap method: {args.bootstrap_method}")

    num_ref_sigs = signatures.shape[1]
    sel_sig_nums = list(range(num_ref_sigs))
    print(f"Analyzing signatures: {signatures.columns[sel_sig_nums].tolist()}")

    # Bring up the cuML GPU backend once, before any solving. Kept behind --use_gpu
    # so the import (cupy/cuml) is only required in GPU processes.
    if getattr(args, 'use_gpu', False):
        configure_gpu_backend(args.gpu_precision)
        print(f"GPU backend enabled (cuML), precision: {args.gpu_precision}")
        print("GPU provenance:", json.dumps(_GPU_SOLVER.provenance(), default=str))

    if args.n_bootstrap > 0:
        # ===== In-process bootstrap loop (P3) =====
        # All bootstrap iterations run in a single process, so Python startup,
        # numba JIT compilation and CSV loading are paid once instead of once per
        # iteration. This also keeps a single (future) GPU context alive across all
        # iterations rather than recreating one per Nextflow task. Each iteration
        # resamples the base table afresh and writes its own indexed output set,
        # for every threshold pair.
        method = args.bootstrap_method

        # Limit samples once; resample the same base table each iteration.
        if args.number_of_samples != -1:
            input_mutations = input_mutations.iloc[:, :args.number_of_samples]

        central = None
        if method == "bootstrap_residuals":
            # the central fit whose residuals every iteration resamples
            central, _ = attribute_at_thresholds(input_mutations.reindex(signatures.index), signatures,
                                                 sel_sig_nums, args, weak_thresholds, strong_thresholds)

        for weak in weak_thresholds:
            for strong in strong_thresholds:
                make_folder_if_not_exists(f'{output_dir(weak, strong)}/bootstrap_output')

        start_index = args.bootstrap_start_index
        end_index = start_index + args.n_bootstrap - 1
        print(f"Running {args.n_bootstrap} bootstrap iteration(s) "
              f"(indices {start_index}..{end_index}, method: {method})")
        start_time = time.process_time()
        for i in range(start_index, end_index + 1):
            if args.seed is not None:
                # Seed per iteration index (not once per process) so iteration i draws the
                # same resample whether iterations are fanned out one-per-job (CPU) or run
                # together in one process (GPU). SeedSequence mixes (seed, i) to avoid the
                # correlation of consecutive integer seeds; it seeds the legacy global RNG
                # that bootstrap_mutation_table() uses.
                iter_seed = np.random.SeedSequence([args.seed, i]).generate_state(1)[0]
                np.random.seed(int(iter_seed))
            results = bootstrap_attributions(input_mutations, signatures, sel_sig_nums, args,
                                             weak_thresholds, strong_thresholds, central)
            for (weak, strong), (b_weights, b_mutations, b_stat_info, _, _) in results.items():
                prefix = (f'{output_dir(weak, strong)}/bootstrap_output/'
                          f'output_{dataset_name}_{mutation_type}{bootstrap_suffix(weak, strong)}_{i}')
                b_weights.to_csv(f'{prefix}_weights_table.csv')
                b_mutations.to_csv(f'{prefix}_mutations_table.csv')
                b_stat_info.to_csv(f'{prefix}_stat_info.csv')
        elapsed = time.process_time() - start_time
        print(f"{args.n_bootstrap} bootstrap iterations took {elapsed:.2f}s "
              f"({elapsed/args.n_bootstrap:.3f}s per iteration)")
    else:
        # ===== Single attribution (central run, or a single externally-indexed bootstrap) =====
        # Limit samples if specified (before resampling, as in the in-process loop)
        if args.number_of_samples != -1:
            input_mutations = input_mutations.iloc[:, :args.number_of_samples]

        # Process samples
        start_time = time.process_time()

        if args.bootstrap:
            # the same resample as iteration `bootstrap_start_index` of the in-process loop
            central = None
            if args.bootstrap_method == "bootstrap_residuals":
                central, _ = attribute_at_thresholds(input_mutations.reindex(signatures.index), signatures,
                                                     sel_sig_nums, args, weak_thresholds, strong_thresholds)
            if args.seed is not None:
                iter_seed = np.random.SeedSequence([args.seed, args.bootstrap_start_index]).generate_state(1)[0]
                np.random.seed(int(iter_seed))
            results = bootstrap_attributions(input_mutations, signatures, sel_sig_nums, args,
                                             weak_thresholds, strong_thresholds, central)
        else:
            input_mutations = input_mutations.reindex(signatures.index)
            trajectories = None
            if auto:
                # choose the weak thresholds from this data's removal trajectories, then reuse them
                trajectories = removal_trajectories(input_mutations, signatures, sel_sig_nums, args)
                weak_thresholds = auto_weak_thresholds([t for _, t in trajectories.values()], args.auto_thresholds)
                print(f"Automatic weak thresholds ({len(trajectories)} samples): "
                      f"{[format_threshold(w) for w in weak_thresholds]}")
                if args.resolve_thresholds_only:
                    write_threshold_labels(args.thresholds_dir, weak_thresholds, strong_thresholds)
                    sys.exit(0)
            results, _ = attribute_at_thresholds(input_mutations, signatures, sel_sig_nums, args,
                                                 weak_thresholds, strong_thresholds, trajectories)
        if multi and args.thresholds_dir:
            write_threshold_labels(args.thresholds_dir, weak_thresholds, strong_thresholds)

        end_time = time.process_time()
        elapsed = end_time - start_time
        n_samples = len(input_mutations.columns)
        print(f"Attribution took {elapsed:.2f}s ({elapsed/n_samples:.2f}s per sample)")

        # Save outputs
        for (weak, strong), outputs in results.items():
            output_path = output_dir(weak, strong)
            make_folder_if_not_exists(output_path)
            output_weights, output_mutations, output_stat_info, residuals_df, fitted_df = outputs
            output_weights.to_csv(f'{output_path}/output_{dataset_name}_{mutation_type}_weights_table.csv')
            output_mutations.to_csv(f'{output_path}/output_{dataset_name}_{mutation_type}_mutations_table.csv')
            output_stat_info.to_csv(f'{output_path}/output_{dataset_name}_{mutation_type}_stat_info.csv')
            residuals_df.to_csv(f'{output_path}/output_{dataset_name}_{mutation_type}_residuals.csv')
            fitted_df.to_csv(f'{output_path}/output_{dataset_name}_{mutation_type}_fitted_values.csv')

            if _GPU_SOLVER is not None:
                # Record the GPU environment next to the attribution it produced
                with open(f'{output_path}/gpu_provenance.json', 'w') as f:
                    json.dump(_GPU_SOLVER.provenance(), f, indent=2, default=str)
