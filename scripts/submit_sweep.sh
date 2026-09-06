#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PARAMS_FILE="${1:-${REPO_ROOT}/parameters/sweep_40.tsv}"
MAX_CONCURRENT="${MAX_CONCURRENT:-8}"

if [[ -z "${RUN_ROOT:-}" && -z "${SCRATCH:-}" ]]; then
  echo "SCRATCH is not set. Set SCRATCH or provide RUN_ROOT explicitly." >&2
  exit 1
fi
if [[ ! -f "${PARAMS_FILE}" ]]; then
  echo "Parameter file not found: ${PARAMS_FILE}" >&2
  exit 1
fi

num_jobs="$(awk '
  /^[[:space:]]*$/ { next }
  /^[[:space:]]*#/ { next }
  { count++ }
  END { print count + 0 }
' "${PARAMS_FILE}")"
if [[ "${num_jobs}" -le 0 ]]; then
  echo "No parameter rows found in ${PARAMS_FILE}" >&2
  exit 1
fi

if [[ -n "${RUN_ROOT:-}" ]]; then
  run_root="${RUN_ROOT}"
else
  mkdir -p "${SCRATCH}/GPES"
  # Atomically reserve a unique directory even for submissions in the same second.
  run_root="$(mktemp -d "${SCRATCH}/GPES/sweep_$(date +%Y%m%d_%H%M%S)_XXXXXX")"
fi
mkdir -p "${run_root}/slurm"
echo "Results:   ${run_root}"

build_job="$(sbatch --parsable \
  --export="ALL,REPO_ROOT=${REPO_ROOT},RUN_ROOT=${run_root}" \
  --output="${run_root}/slurm/build_%j.out" \
  --error="${run_root}/slurm/build_%j.err" \
  "${SCRIPT_DIR}/build.slurm")"
array_job="$(sbatch --parsable \
  --export="ALL,REPO_ROOT=${REPO_ROOT},RUN_ROOT=${run_root}" \
  --dependency="afterok:${build_job}" \
  --array="1-${num_jobs}%${MAX_CONCURRENT}" \
  --output="${run_root}/slurm/%A_%a.out" \
  --error="${run_root}/slurm/%A_%a.err" \
  "${SCRIPT_DIR}/run.slurm" "${PARAMS_FILE}")"

echo "Build job: ${build_job}"
echo "Sweep job: ${array_job} (${num_jobs} tasks, at most ${MAX_CONCURRENT} concurrent)"
