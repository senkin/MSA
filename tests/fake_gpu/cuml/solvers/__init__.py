"""SciPy stand-in for cuml.solvers.nnls_batched (rapidsai/cuml PR #8402).

Same call signature and layout as the real solver: A shared, B gathered by b_index,
column-major masks (n_signatures, n_problems). Each problem is solved with scipy in
fp64 and the result is returned in the input dtype, so 'mixed'/'fp32' runs see fp32
storage rounding (not an fp32 solve). This checks MSA's wrapper, not cuML itself.
"""
import numpy as np
from scipy.optimize import nnls


def nnls_batched(A, B, masks=None, b_index=None, *, solver='lawson', compute_fitted=True,
                 maxiter=None, tol=1e-6):
    # signature of the pinned cuML commit (see CUML_COMMIT in setup_gpu_env.sh)
    if masks is None:
        masks = np.ones((A.shape[1], B.shape[1]), dtype=bool)
    n_problems = masks.shape[1]
    if b_index is None:
        b_index = np.arange(n_problems)
    A64 = np.asarray(A, dtype=np.float64)
    B64 = np.asarray(B, dtype=np.float64)
    X = np.zeros((A.shape[1], n_problems))
    fitted = np.zeros((A.shape[0], n_problems))
    for j in range(n_problems):
        cols = np.flatnonzero(masks[:, j])
        w, _ = nnls(A64[:, cols], B64[:, int(b_index[j])])
        X[cols, j] = w
        fitted[:, j] = A64[:, cols] @ w
    X = X.astype(A.dtype, copy=False)
    fitted = fitted.astype(A.dtype, copy=False)
    return (X, fitted) if compute_fitted else (X, None)
