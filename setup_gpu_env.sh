#!/usr/bin/env bash
#
# setup_gpu_env.sh - one-time provisioning of the GPU (cuML) conda env for MSA.
#
# Builds the batched-NNLS cuML (rapidsai/cuml PR #8402, not yet released) from
# source into a conda env, then tops it up with MSA's runtime deps. Run this ONCE
# per cluster + GPU architecture, as a build job on a node that has an NVIDIA driver
# and enough CPU/RAM - the CUDA/C++ compile is heavy (tens of GB RAM, ~30+ min).
# Do NOT run it on a login/head node. The CUDA toolkit (nvcc) comes from the conda
# env, so no system CUDA install is needed, but a matching NVIDIA driver must exist.
#
# When it finishes, set in nextflow.config (or pass --gpu_conda_env) to the FULL
# ENV PATH it prints (e.g. ~/miniforge3/envs/all_cuda-133_arch-x86_64) - a bare name
# would be treated by Nextflow's conda directive as a package to install.
#
# It is idempotent: an existing clone or conda env is reused (delete them to rebuild).

set -euo pipefail

usage() {
  sed -n '3,15p' "$0" | sed 's/^# \{0,1\}//'
  cat <<'EOF'

Usage: ./setup_gpu_env.sh [options]
  -d, --build-dir DIR   where to clone cuml         (default: $PWD/cuml-build)
  -e, --env-name NAME   conda env name to create    (default: <cuda-env>)
  -c, --cuda-env NAME   RAPIDS env file basename     (default: all_cuda-133_arch-x86_64)
  -b, --branch NAME     cuml branch                  (default: fea-nnls)
  -r, --repo URL        cuml repo                    (default: achirkin/cuml)
      --allgpuarch      build for all GPU archs (portable, slower) instead of native
      --no-overlay      skip applying environment_gpu.yml
  -h, --help            show this help
EOF
}

# --- defaults (override via flags or env) -----------------------------------
REPO_URL="${CUML_REPO_URL:-https://github.com/achirkin/cuml.git}"
BRANCH="${CUML_BRANCH:-fea-nnls}"
CUDA_ENV="${CUML_CUDA_ENV:-all_cuda-133_arch-x86_64}"   # conda/environments/<CUDA_ENV>.yaml
ENV_NAME="${GPU_ENV_NAME:-}"                            # defaults to $CUDA_ENV below
BUILD_DIR="${CUML_BUILD_DIR:-$PWD/cuml-build}"
ALLGPUARCH=0
APPLY_OVERLAY=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OVERLAY_YAML="$SCRIPT_DIR/environment_gpu.yml"

# --- parse flags ------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -d|--build-dir) BUILD_DIR="$2"; shift 2 ;;
    -e|--env-name)  ENV_NAME="$2"; shift 2 ;;
    -c|--cuda-env)  CUDA_ENV="$2"; shift 2 ;;
    -b|--branch)    BRANCH="$2"; shift 2 ;;
    -r|--repo)      REPO_URL="$2"; shift 2 ;;
    --allgpuarch)   ALLGPUARCH=1; shift ;;
    --no-overlay)   APPLY_OVERLAY=0; shift ;;
    -h|--help)      usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done
ENV_NAME="${ENV_NAME:-$CUDA_ENV}"

# --- prerequisites ----------------------------------------------------------
command -v git >/dev/null || { echo "error: git not found" >&2; exit 1; }
CONDA_EXE="${CONDA_EXE:-conda}"
command -v "$CONDA_EXE" >/dev/null || { echo "error: conda not found (activate a base conda first)" >&2; exit 1; }
# RAPIDS envs solve painfully slowly with classic conda; prefer mamba when present.
SOLVER_EXE="$CONDA_EXE"
command -v mamba >/dev/null && SOLVER_EXE="mamba"
# make 'conda activate' work in this non-interactive shell
CONDA_BASE="$("$CONDA_EXE" info --base)"
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"

cat <<EOF
>>> GPU (cuML) env setup for MSA
    repo      : $REPO_URL @ $BRANCH
    cuda env  : conda/environments/$CUDA_ENV.yaml
    env name  : $ENV_NAME
    build dir : $BUILD_DIR
    gpu archs : $([[ $ALLGPUARCH == 1 ]] && echo 'all (portable)' || echo 'native (this node)')
    solver    : $SOLVER_EXE
    overlay   : $([[ $APPLY_OVERLAY == 1 ]] && echo "$OVERLAY_YAML" || echo skip)
EOF

# --- 1. clone (or reuse) the branch ----------------------------------------
if [[ -d "$BUILD_DIR/.git" ]]; then
  echo ">>> Reusing existing clone at $BUILD_DIR"
else
  echo ">>> Cloning $REPO_URL ($BRANCH)"
  git clone --depth 1 -b "$BRANCH" "$REPO_URL" "$BUILD_DIR"
fi
cd "$BUILD_DIR"

ENV_YAML="conda/environments/$CUDA_ENV.yaml"
[[ -f "$ENV_YAML" ]] || { echo "error: env file not found: $BUILD_DIR/$ENV_YAML" >&2; exit 1; }

# --- 2. create (or reuse) the conda env ------------------------------------
if "$CONDA_EXE" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo ">>> Conda env '$ENV_NAME' already exists, reusing (delete it to rebuild from scratch)"
else
  echo ">>> Creating conda env '$ENV_NAME' from $ENV_YAML"
  "$SOLVER_EXE" env create -n "$ENV_NAME" -f "$ENV_YAML"
fi

# conda's activation hooks are not 'set -u' clean - the cuda-nvcc hook references
# NVCC_PREPEND_FLAGS without a default, which aborts under nounset. Relax it here.
set +u
conda activate "$ENV_NAME"
set -u
# Full env path - this, NOT the bare name, is what gpu_conda_env must be set to:
# Nextflow's conda directive treats a bare name as a package to install.
ENV_PREFIX="${CONDA_PREFIX:-$ENV_NAME}"

# --- 3. build + install cuML into the env ----------------------------------
echo ">>> Building libcuml + cuml (this is the slow part)"
if [[ $ALLGPUARCH == 1 ]]; then
  ./build.sh --allgpuarch libcuml cuml
else
  ./build.sh libcuml cuml
fi

# --- 4. optional MSA overlay -----------------------------------------------
if [[ $APPLY_OVERLAY == 1 && -f "$OVERLAY_YAML" ]]; then
  echo ">>> Applying MSA overlay ($OVERLAY_YAML)"
  "$SOLVER_EXE" env update -n "$ENV_NAME" -f "$OVERLAY_YAML"
fi

# --- 5. verify --------------------------------------------------------------
echo ">>> Verifying imports"
python - <<'PY'
import importlib, sys
ok = True
# MSA host deps that run_NNLS.py / common_methods.py import at load - must be present.
for m in ["numpy", "pandas", "scipy", "numba",
          "matplotlib.pyplot", "statsmodels.stats.proportion"]:
    try:
        importlib.import_module(m)
    except Exception as e:
        ok = False
        print(f"  MISSING MSA dep {m}: {e}", file=sys.stderr)
# GPU stack needs a driver/GPU at import; a failure here on a CPU-only build node is
# expected - re-run the check on a GPU node before using the pipeline.
gpu_ok = True
for m in ["cupy", "rmm", "cuml"]:
    try:
        importlib.import_module(m)
    except Exception as e:
        gpu_ok = False
        print(f"  NOTE: could not import {m} here ({type(e).__name__}); "
              "verify on a GPU node", file=sys.stderr)
if gpu_ok:
    from cuml.solvers import nnls_batched
    print("  cuml.solvers.nnls_batched present:", callable(nnls_batched))
if not ok:
    sys.exit(1)
print("MSA deps OK." + ("" if gpu_ok else " (GPU stack unverified on this node)"))
PY

cat <<EOF

>>> Done. In nextflow.config (or via --gpu_conda_env) set:
        use_GPU       = true
        gpu_conda_env = '$ENV_PREFIX'
    (use this full path, not the bare env name - Nextflow's conda directive treats
    a bare name as a package to install). Run with a conda-enabled profile
    (e.g. -profile conda) so the labelled NNLS processes pick up this env.
EOF
