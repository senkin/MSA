import argparse
import os
import warnings
import pandas as pd
import numpy as np
import time
from scipy.optimize import nnls
from numba import jit
from common_methods import make_folder_if_not_exists, calculate_similarity, format_threshold

# Vectorized bootstrap methods
def bootstrap_mutation_table(input_dataframe, method="classic", fitted=None, residuals=None):
    """Optimized bootstrap with vectorized operations."""
    if method == "classic":
        return input_dataframe.sample(n=len(input_dataframe), replace=True)

    elif method == "bootstrap_residuals":
        if residuals is None or fitted is None:
            raise ValueError(f"Residuals and fitted values required for {method}")
        return fitted + bootstrap_mutation_table(residuals, method="classic")

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
    bytes of signal (a 4608:1 amplification at SBS-4608), so a GPU backend should
    fuse the reduction and implement batched_solve_and_score() instead.

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
    reduction and keep both on the device, returning only ``sims``. That is worth
    a great deal at high context: handing back ``fitted`` instead would move
    n_channels * 8 bytes per problem (~1.8 TB over a 1000-bootstrap SBS-4608 run)
    to deliver 8 bytes of signal per problem (~0.4 GB).

    ``A``, ``B`` and ``norm_obs`` are loop-invariant across the greedy steps of a
    given batch, so a backend is free to keep them device-resident between calls.

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


def process_samples_batch(input_mutations, signatures, sel_sig_nums, args):
    """Process all samples, accumulating results in preallocated numpy arrays.

    Per-sample results are written into numpy arrays and assembled into output
    DataFrames once at the end. The original implementation assigned into
    DataFrames row-by-row via .loc, which is slow at scale.
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

            per_sample_cols = optimise_signatures_batched(
                B, norm_obs, sig_values, initial_cols,
                strategy=args.optimisation_strategy,
                weak_threshold=args.weak_threshold,
                chunk_size=getattr(args, 'gpu_batch_size', 4096))

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

            if args.optimise_signatures:
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
    parser.add_argument("-W", "--weak_threshold", dest="weak_threshold", default=0.01, type=float)
    parser.add_argument("-S", "--strong_threshold", dest="strong_threshold", default=0.01, type=float)
    parser.add_argument("-N", "--normalise_mutations", dest="normalise_mutations", action="store_true")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true")
    parser.add_argument("-n", "--number", dest="number_of_samples", default=-1, type=int)
    parser.add_argument("-B", "--bootstrap", dest="bootstrap", action="store_true")
    parser.add_argument("--bootstrap_method", dest="bootstrap_method", default='binomial')
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
    parser.add_argument("--bootstrap_output_suffix", dest="bootstrap_output_suffix", default='',
                       help="String inserted into per-iteration bootstrap output filenames "
                            "(e.g. '<weak>_<strong>' thresholds), to match downstream expectations.")
    parser.add_argument("--use_gpu", dest="use_gpu", action="store_true",
                       help="Use the sample-transposed batched optimisation path (removal strategy), "
                            "which issues one large masked NNLS batch per greedy step across all "
                            "still-active samples. Intended for GPU offloading via "
                            "batched_solve_and_score().")
    parser.add_argument("--gpu_batch_size", dest="gpu_batch_size", default=4096, type=int,
                       help="Maximum number of NNLS problems solved per batched call (memory bound).")

    args = parser.parse_args()

    # Setup paths
    dataset_name = args.dataset_name
    mutation_type = args.mutation_type
    context = args.context
    output_path = f"{args.output_path}/{dataset_name}"

    if args.add_suffix:
        output_path += f'_{context}_NNLS'
        if args.optimise_signatures:
            weak_str = format_threshold(args.weak_threshold)
            strong_str = format_threshold(args.strong_threshold)
            output_path += f'_{weak_str}_{strong_str}'
        else:
            output_path += '_unoptimised'

    make_folder_if_not_exists(output_path)

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
        print(f"Optimised NNLS: weak/strong thresholds = {args.weak_threshold}/{args.strong_threshold}")
    if args.bootstrap:
        print(f"Bootstrap method: {args.bootstrap_method}")

    num_ref_sigs = signatures.shape[1]
    sel_sig_nums = list(range(num_ref_sigs))
    print(f"Analyzing signatures: {signatures.columns[sel_sig_nums].tolist()}")

    if args.n_bootstrap > 0:
        # ===== In-process bootstrap loop (P3) =====
        # All bootstrap iterations run in a single process, so Python startup,
        # numba JIT compilation and CSV loading are paid once instead of once per
        # iteration. This also keeps a single (future) GPU context alive across all
        # iterations rather than recreating one per Nextflow task. Each iteration
        # resamples the base table afresh and writes its own indexed output set.
        method = args.bootstrap_method

        # Limit samples once; resample the same base table each iteration.
        if args.number_of_samples != -1:
            input_mutations = input_mutations.iloc[:, :args.number_of_samples]

        residuals_table = fitted_table = None
        if method == "bootstrap_residuals":
            residuals_table = pd.read_csv(
                f'{inp_path}/{dataset_name}/output_{dataset_name}_{mutation_type}_residuals.csv',
                index_col=[0, 1, 2])
            fitted_table = pd.read_csv(
                f'{inp_path}/{dataset_name}/output_{dataset_name}_{mutation_type}_fitted_values.csv',
                index_col=[0, 1, 2])

        bootstrap_dir = f'{output_path}/bootstrap_output'
        make_folder_if_not_exists(bootstrap_dir)
        suffix = f'_{args.bootstrap_output_suffix}' if args.bootstrap_output_suffix else ''

        start_index = args.bootstrap_start_index
        end_index = start_index + args.n_bootstrap - 1
        print(f"Running {args.n_bootstrap} bootstrap iteration(s) "
              f"(indices {start_index}..{end_index}, method: {method})")
        start_time = time.process_time()
        for i in range(start_index, end_index + 1):
            if method == "bootstrap_residuals":
                resampled = bootstrap_mutation_table(input_mutations, method=method,
                                                     fitted=fitted_table, residuals=residuals_table)
            else:
                resampled = bootstrap_mutation_table(input_mutations, method=method)
            resampled = resampled.reindex(signatures.index)

            b_weights, b_mutations, b_stat_info, _, _ = \
                process_samples_batch(resampled, signatures, sel_sig_nums, args)

            prefix = f'{bootstrap_dir}/output_{dataset_name}_{mutation_type}{suffix}_{i}'
            b_weights.to_csv(f'{prefix}_weights_table.csv')
            b_mutations.to_csv(f'{prefix}_mutations_table.csv')
            b_stat_info.to_csv(f'{prefix}_stat_info.csv')
        elapsed = time.process_time() - start_time
        print(f"{args.n_bootstrap} bootstrap iterations took {elapsed:.2f}s "
              f"({elapsed/args.n_bootstrap:.3f}s per iteration)")
    else:
        # ===== Single attribution (central run, or a single externally-indexed bootstrap) =====
        # Bootstrap if requested
        if args.bootstrap and args.bootstrap_method != "bootstrap_residuals":
            input_mutations = bootstrap_mutation_table(input_mutations, method=args.bootstrap_method)

        # Limit samples if specified
        if args.number_of_samples != -1:
            input_mutations = input_mutations.iloc[:, :args.number_of_samples]

        # Align indices
        input_mutations = input_mutations.reindex(signatures.index)

        # Handle bootstrap_residuals method
        if args.bootstrap and args.bootstrap_method == "bootstrap_residuals":
            residuals_df = pd.read_csv(f'{inp_path}/{dataset_name}/output_{dataset_name}_{mutation_type}_residuals.csv',
                                       index_col=[0,1,2])
            fitted_df = pd.read_csv(f'{inp_path}/{dataset_name}/output_{dataset_name}_{mutation_type}_fitted_values.csv',
                                   index_col=[0,1,2])
            input_mutations = bootstrap_mutation_table(input_mutations, method=args.bootstrap_method,
                                                       fitted=fitted_df, residuals=residuals_df)

        # Process samples
        start_time = time.process_time()

        output_weights, output_mutations, output_stat_info, residuals_df, fitted_df = \
            process_samples_batch(input_mutations, signatures, sel_sig_nums, args)

        end_time = time.process_time()
        elapsed = end_time - start_time
        n_samples = len(input_mutations.columns)
        print(f"Attribution took {elapsed:.2f}s ({elapsed/n_samples:.2f}s per sample)")

        # Save outputs
        output_weights.to_csv(f'{output_path}/output_{dataset_name}_{mutation_type}_weights_table.csv')
        output_mutations.to_csv(f'{output_path}/output_{dataset_name}_{mutation_type}_mutations_table.csv')
        output_stat_info.to_csv(f'{output_path}/output_{dataset_name}_{mutation_type}_stat_info.csv')
        residuals_df.to_csv(f'{output_path}/output_{dataset_name}_{mutation_type}_residuals.csv')
        fitted_df.to_csv(f'{output_path}/output_{dataset_name}_{mutation_type}_fitted_values.csv')

        if not args.bootstrap:
            residuals_df.to_csv(f'{inp_path}/{dataset_name}/output_{dataset_name}_{mutation_type}_residuals.csv')
            fitted_df.to_csv(f'{inp_path}/{dataset_name}/output_{dataset_name}_{mutation_type}_fitted_values.csv')
