# GPU backend state. _GPU_CACHE keeps A / B / norm_obs resident on the device
# for one "generation" (one (A, B, norm_obs) identity triple); they are
# loop-invariant across the greedy steps of an optimisation batch, so this
# uploads each at most once per run and reuses it across calls.
_GPU_INITIALISED = False
_GPU_CACHE = {}
# Kept alive for the process lifetime: the per-thread default stream that both
# cupy and the cuML/raft handle share (see _init_gpu_backend).
_GPU_STREAM = None


def _init_gpu_backend():
    """Enable a bounded RMM pool and route cupy through it (once).

    The greedy optimisation issues many short-lived device allocations; the pool
    avoids the per-call cudaMalloc/cudaFree churn that would otherwise dominate.

    Pool size is set by the launcher via ``RMM_POOL_SIZE_MIB`` (MiB). When unset,
    default to 90% of total GPU memory (single-process). The pool grows on demand
    up to that cap — it is not pre-allocated.

    Also switches cupy onto the per-thread default stream so it shares the exact
    stream the cuML/raft handle uses (raft's ``Handle()`` is built on
    ``cudaStreamPerThread``). Otherwise cupy runs on the legacy default stream,
    which does not synchronise with the raft stream: the Lawson kernel could then
    read the H2D copies / gathers cupy is still producing, an unordered data race
    that surfaces intermittently as CUDA_ERROR_ILLEGAL_ADDRESS at handle.sync().
    """
    global _GPU_INITIALISED, _GPU_STREAM
    if _GPU_INITIALISED:
        return
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator
    import cupy as cp
    # Pool cap comes from the launcher (MiB); it already accounts for the number
    # of concurrent processes sharing this GPU. Fall back to 60% of total memory
    # (single-process) when unset. Grow on demand — do not pre-allocate.
    env_mib = os.environ.get("RMM_POOL_SIZE_MIB")
    if env_mib is not None:
        pool_size = int(float(env_mib) * (1 << 20))
    else:
        _, total_memory = cp.cuda.runtime.memGetInfo()
        pool_size = int(total_memory * 0.6)
    rmm.reinitialize(
        pool_allocator=True,
        initial_pool_size=0,
        maximum_pool_size=pool_size,
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)
    # Match the raft handle's stream (cudaStreamPerThread) for all cupy ops so
    # inputs to the cuML kernels are ordered on the same stream. Held in a global
    # so the active stream is not garbage-collected mid-run.
    _GPU_STREAM = cp.cuda.Stream(ptds=True)
    _GPU_STREAM.use()
    _GPU_INITIALISED = True


def _gpu_resident(A, B, norm_obs):
    """Return device (cupy) copies of the loop-invariant A, B, norm_obs.

    Cached by the identity triple so repeated calls in one optimisation batch
    reuse the same device arrays; a new triple evicts the previous generation to
    keep device memory bounded.
    """
    import cupy as cp
    gen = (id(A), id(B), id(norm_obs))
    if _GPU_CACHE.get("gen") != gen:
        _GPU_CACHE.clear()
        _GPU_CACHE["gen"] = gen
        # Pin the host objects so their id() cannot be reused while cached.
        _GPU_CACHE["refs"] = (A, B, norm_obs)
        _GPU_CACHE["A"] = cp.asfortranarray(cp.asarray(A))
        _GPU_CACHE["B"] = cp.asfortranarray(cp.asarray(B))
        _GPU_CACHE["norm"] = cp.asarray(norm_obs)
    return _GPU_CACHE["A"], _GPU_CACHE["B"], _GPU_CACHE["norm"]


def _gpu_batched_solve_and_score(A, B, norm_obs, masks, b_index, chunk_size):
    """GPU implementation of batched_solve_and_score().

    Solves each masked NNLS problem with cuml.solvers.nnls_batched and reduces
    to the L2_normalised_by_first similarity entirely on the device, so only the
    (n_problems,) similarity vector is copied back to the host. ``fitted`` never
    leaves the GPU.
    """
    import cupy as cp
    from cuml.solvers import nnls_batched as _cuml_nnls_batched

    A_dev, B_dev, norm_dev = _gpu_resident(A, B, norm_obs)

    n_problems = masks.shape[1]
    sims = np.empty(n_problems)
    for start in range(0, n_problems, chunk_size):
        stop = min(start + chunk_size, n_problems)
        # masks is F-contiguous, so a column block is itself F-contiguous and
        # rides to the device (and through cuML) without a layout copy.
        masks_dev = cp.asarray(masks[:, start:stop])
        bi_dev = cp.asarray(b_index[start:stop])
        _, fitted = _cuml_nnls_batched(
            A_dev, B_dev, masks_dev, b_index=bi_dev,
            solver=_NNLS_SOLVER, compute_fitted=True)
        resid = cp.linalg.norm(B_dev[:, bi_dev] - fitted, axis=0)
        sims[start:stop] = cp.asnumpy(1.0 - resid / norm_dev[bi_dev])
    return sims
