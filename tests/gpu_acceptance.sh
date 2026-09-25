#!/usr/bin/env bash
#
# gpu_acceptance.sh - accept (or reject) a GPU/cuML environment for MSA.
#
# Attributes the first --samples samples of SIM_test_GPU (SBS96, 65 signatures,
# optimised removal at --threshold) plus --bootstraps seeded bootstrap iterations, on
# the CPU (fp64 reference) and on the GPU in several configurations, then compares the
# attributions with bin/compare_attributions.py:
#
#   check                              expectation                        gating
#   GPU fp64        vs CPU             <= --max-differing-fraction        yes
#   GPU fp64 repeat vs GPU fp64        identical (determinism)            yes
#   GPU fp64 batch  vs GPU fp64        identical (--gpu_batch_size 4096   yes
#                                      vs auto: results must not depend
#                                      on the GPU's memory size)
#   GPU mixed       vs GPU fp64        <= --max-differing-fraction        yes (default precision)
#   GPU fp32        vs GPU fp64        reported                           no
#   GPU trajectory  vs one run per     identical (the default multi-      yes
#   (all thresholds threshold          threshold scan replays one batched
#   in one run)                        trajectory; requires per-problem
#                                      results independent of batching)
#
# Samples whose greedy decision sits within rounding error of the threshold can flip
# between precisions or solvers, hence a small allowed fraction rather than equality.
#
# Run on a GPU node, e.g.:
#   tests/gpu_acceptance.sh --gpu-python ~/miniforge3/envs/all_cuda-133_arch-x86_64/bin/python
# --fake-gpu runs the GPU configurations on the CPU through tests/fake_gpu (checks this
# script, not the GPU).

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GPU_PYTHON=""
CPU_PYTHON="python"
SAMPLES=200
BOOTSTRAPS=4
THRESHOLD=0.0001
SEED=1
WORKDIR="$PWD/gpu_acceptance"
CPU_JOBS="$(nproc 2>/dev/null || echo 4)"
TOLERANCE=0.001
MAX_FRACTION=0.02
FAKE_GPU=0

usage() {
  sed -n '3,25p' "$0" | sed 's/^# \{0,1\}//'
  cat <<EOF

Usage: tests/gpu_acceptance.sh --gpu-python PATH [options]
  --gpu-python PATH            python of the GPU (cuML) env                      (required unless --fake-gpu)
  --cpu-python PATH            python of the CPU MSA env                         (default: $CPU_PYTHON)
  --samples N                  samples attributed                                (default: $SAMPLES)
  --bootstraps K               seeded bootstrap iterations                       (default: $BOOTSTRAPS)
  --threshold W                weak threshold                                    (default: $THRESHOLD)
  --seed S                     bootstrap seed                                    (default: $SEED)
  --workdir DIR                outputs and logs                                  (default: $WORKDIR)
  --cpu-jobs J                 parallel CPU processes for the reference          (default: $CPU_JOBS)
  --tolerance T                weight difference counted as a mismatch           (default: $TOLERANCE)
  --max-differing-fraction F   allowed fraction of mismatching attributions      (default: $MAX_FRACTION)
  --fake-gpu                   run the GPU configurations on tests/fake_gpu (CPU)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gpu-python) GPU_PYTHON="$2"; shift 2 ;;
    --cpu-python) CPU_PYTHON="$2"; shift 2 ;;
    --samples) SAMPLES="$2"; shift 2 ;;
    --bootstraps) BOOTSTRAPS="$2"; shift 2 ;;
    --threshold) THRESHOLD="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --workdir) WORKDIR="$2"; shift 2 ;;
    --cpu-jobs) CPU_JOBS="$2"; shift 2 ;;
    --tolerance) TOLERANCE="$2"; shift 2 ;;
    --max-differing-fraction) MAX_FRACTION="$2"; shift 2 ;;
    --fake-gpu) FAKE_GPU=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

GPU_ENV=()
if [[ $FAKE_GPU == 1 ]]; then
  GPU_PYTHON="$CPU_PYTHON"
  GPU_ENV=(env "PYTHONPATH=$REPO/tests/fake_gpu${PYTHONPATH:+:$PYTHONPATH}")
elif [[ -z "$GPU_PYTHON" ]]; then
  echo "error: --gpu-python is required (or --fake-gpu)" >&2; usage >&2; exit 2
fi

DATASET=SIM_test_GPU
mkdir -p "$WORKDIR/inputs/$DATASET" "$WORKDIR/logs"
WORKDIR="$(cd "$WORKDIR" && pwd)"
# run_NNLS.py writes residuals next to its inputs, so work on a copy
cp "$REPO/input_mutation_tables/$DATASET/WGS_$DATASET.96.csv" "$WORKDIR/inputs/$DATASET/"
export MPLCONFIGDIR="$WORKDIR/.matplotlib"

COMMON=(-d "$DATASET" -t SBS -c 96 -p sigProfiler -i "$WORKDIR/inputs" -s "$REPO/signature_tables"
        -x -W "$THRESHOLD" -S 0 -n "$SAMPLES")
BOOT=(-B --bootstrap_method binomial --seed "$SEED")

run() {  # run <log name> <command...>; show the log tail on failure
  local log="$WORKDIR/logs/$1.log"; shift
  if ! "$@" > "$log" 2>&1; then
    echo "FAILED: $* (log: $log)" >&2; tail -20 "$log" >&2; exit 2
  fi
}

elapsed() { echo $(( $(date +%s) - $1 ))s; }

echo ">>> CPU reference: $SAMPLES samples + $BOOTSTRAPS bootstrap iterations ($CPU_JOBS parallel jobs)"
t0=$(date +%s)
rm -rf "$WORKDIR/cpu"
run cpu_central "$CPU_PYTHON" "$REPO/bin/run_NNLS.py" "${COMMON[@]}" -o "$WORKDIR/cpu" &
central_pid=$!
# iteration i is seeded from (seed, i), so one process per iteration reproduces a single run
seq 1 "$BOOTSTRAPS" | xargs -P "$CPU_JOBS" -I{} "$CPU_PYTHON" "$REPO/bin/run_NNLS.py" "${COMMON[@]}" "${BOOT[@]}" \
    --n_bootstrap 1 --bootstrap_start_index {} -o "$WORKDIR/cpu" > "$WORKDIR/logs/cpu_bootstrap.log" 2>&1 \
  || { echo "FAILED: CPU bootstrap (log: $WORKDIR/logs/cpu_bootstrap.log)" >&2; tail -20 "$WORKDIR/logs/cpu_bootstrap.log" >&2; exit 2; }
wait "$central_pid" || exit 2
echo "    done in $(elapsed "$t0")"

gpu_run() {  # gpu_run <name> <precision> <batch size>
  local name=$1 precision=$2 batch=$3 t0
  t0=$(date +%s)
  rm -rf "$WORKDIR/$name"
  local gpu=(--use_gpu --gpu_precision "$precision" --gpu_batch_size "$batch")
  run "${name}_central" ${GPU_ENV[@]+"${GPU_ENV[@]}"} "$GPU_PYTHON" "$REPO/bin/run_NNLS.py" "${COMMON[@]}" "${gpu[@]}" -o "$WORKDIR/$name"
  run "${name}_bootstrap" ${GPU_ENV[@]+"${GPU_ENV[@]}"} "$GPU_PYTHON" "$REPO/bin/run_NNLS.py" "${COMMON[@]}" "${BOOT[@]}" \
      --n_bootstrap "$BOOTSTRAPS" --bootstrap_start_index 1 "${gpu[@]}" -o "$WORKDIR/$name"
  echo "    $name: $(elapsed "$t0")"
}

echo ">>> GPU runs$([[ $FAKE_GPU == 1 ]] && echo ' (FAKE GPU: tests/fake_gpu on the CPU)')"
gpu_run gpu_fp64 fp64 0
gpu_run gpu_fp64_repeat fp64 0
gpu_run gpu_fp64_batch4096 fp64 4096
gpu_run gpu_mixed mixed 0
gpu_run gpu_fp32 fp32 0

# Multi-threshold scan (default precision): one run replaying the batched trajectory at
# every threshold vs one run per threshold, same output folder layout (--add_suffix)
SCAN=(0 "$THRESHOLD" 0.001 0.01)
t0=$(date +%s)
rm -rf "$WORKDIR/gpu_scan_trajectory" "$WORKDIR/gpu_scan_separate"
GPU_MIXED=(--use_gpu --gpu_precision mixed --gpu_batch_size 0 --add_suffix -S 0)
run scan_trajectory_central ${GPU_ENV[@]+"${GPU_ENV[@]}"} "$GPU_PYTHON" "$REPO/bin/run_NNLS.py" "${COMMON[@]}" \
    "${GPU_MIXED[@]}" -W "${SCAN[@]}" -o "$WORKDIR/gpu_scan_trajectory"
run scan_trajectory_bootstrap ${GPU_ENV[@]+"${GPU_ENV[@]}"} "$GPU_PYTHON" "$REPO/bin/run_NNLS.py" "${COMMON[@]}" \
    "${BOOT[@]}" --n_bootstrap "$BOOTSTRAPS" "${GPU_MIXED[@]}" -W "${SCAN[@]}" -o "$WORKDIR/gpu_scan_trajectory"
for w in "${SCAN[@]}"; do
  run "scan_separate_central_$w" ${GPU_ENV[@]+"${GPU_ENV[@]}"} "$GPU_PYTHON" "$REPO/bin/run_NNLS.py" "${COMMON[@]}" \
      "${GPU_MIXED[@]}" -W "$w" -o "$WORKDIR/gpu_scan_separate"
  run "scan_separate_bootstrap_$w" ${GPU_ENV[@]+"${GPU_ENV[@]}"} "$GPU_PYTHON" "$REPO/bin/run_NNLS.py" "${COMMON[@]}" \
      "${BOOT[@]}" --n_bootstrap "$BOOTSTRAPS" "${GPU_MIXED[@]}" -W "$w" -o "$WORKDIR/gpu_scan_separate" \
      --bootstrap_output_suffix "$("$CPU_PYTHON" -c "import sys; sys.path.insert(0, '$REPO/bin'); from common_methods import format_threshold as f; print(f(float('$w')) + '_0')")"
done
echo "    threshold scan (trajectory + ${#SCAN[@]} separate runs): $(elapsed "$t0")"

echo ">>> GPU environment (from gpu_fp64/$DATASET/gpu_provenance.json)"
"$CPU_PYTHON" - "$WORKDIR/gpu_fp64/$DATASET/gpu_provenance.json" <<'PY'
import json, sys
p = json.load(open(sys.argv[1]))
dev, cuda, ver = p.get('device', {}), p.get('cuda', {}), p.get('versions', {})
print(f"    device: {dev.get('name')} (cc {dev.get('compute_capability')}), "
      f"driver {cuda.get('driver_version')}, runtime {cuda.get('runtime_version')}")
print(f"    cuml {ver.get('cuml')}, cupy {ver.get('cupy')}, rmm {ver.get('rmm')}; auto batch {p.get('auto_batch_size')}")
build = p.get('cuml_build')
print(f"    cuml build: {build.get('commit') if isinstance(build, dict) else build}")
PY

compare() {  # compare <label> <gating 0/1> <reference> <candidate> [compare_attributions options]
  local label=$1 gating=$2 ref=$3 cand=$4; shift 4
  echo ">>> $label"
  if "$CPU_PYTHON" "$REPO/bin/compare_attributions.py" "$WORKDIR/$ref" "$WORKDIR/$cand" \
        --json "$WORKDIR/logs/compare_${cand}_vs_${ref}.json" "$@" | sed 's/^/    /'; then
    RESULTS+=("PASS    $label")
  elif [[ $gating == 1 ]]; then
    RESULTS+=("FAIL    $label"); FAILED=1
  else
    RESULTS+=("REPORT  $label (outside tolerance, not gating)")
  fi
}

RESULTS=(); FAILED=0
TOL=(--tolerance "$TOLERANCE" --max-differing-fraction "$MAX_FRACTION")
compare "GPU fp64 vs CPU (fp64 reference)" 1 cpu gpu_fp64 "${TOL[@]}"
compare "GPU fp64 repeat vs GPU fp64 (determinism)" 1 gpu_fp64 gpu_fp64_repeat --exact
compare "GPU fp64 batch 4096 vs auto batch (batch-size invariance)" 1 gpu_fp64 gpu_fp64_batch4096 --exact
compare "GPU mixed vs GPU fp64 (default precision)" 1 gpu_fp64 gpu_mixed "${TOL[@]}"
compare "GPU fp32 vs GPU fp64" 0 gpu_fp64 gpu_fp32 "${TOL[@]}"
compare "GPU threshold scan: trajectory vs one run per threshold" 1 gpu_scan_separate gpu_scan_trajectory --exact

echo
echo "================ GPU acceptance summary ================"
printf '%s\n' "${RESULTS[@]}"
echo "Outputs and per-check JSON reports: $WORKDIR"
exit $FAILED
