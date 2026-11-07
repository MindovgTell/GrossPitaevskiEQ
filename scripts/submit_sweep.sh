#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PARAMS_FILE="${1:-${REPO_ROOT}/parameters/sweep_80.tsv}"
MAX_CONCURRENT="${MAX_CONCURRENT:-8}"

if [[ -z "${RUN_ROOT:-}" && -z "${SCRATCH:-}" ]]; then
  echo "SCRATCH is not set. Set SCRATCH or provide RUN_ROOT explicitly." >&2
  exit 1
fi
if [[ ! -f "${PARAMS_FILE}" ]]; then
  echo "Parameter file not found: ${PARAMS_FILE}" >&2
  exit 1
fi
if ! awk '
  /^[[:space:]]*$/ || /^[[:space:]]*#/ { next }
  {
    count++
    if (NF < 8 || NF > 11) {
      print "Parameter row " count " requires 8 to 11 fields" > "/dev/stderr"
      invalid = 1
      next
    }
    name = NF >= 9 ? $9 : sprintf("run_%04d", count)
    if (name !~ /^[A-Za-z0-9][A-Za-z0-9._-]*$/ || seen[name]++) {
      print "Invalid or repeated run name in row " count ": " name > "/dev/stderr"
      invalid = 1
    }
  }
  END { if (invalid || count == 0) exit 1 }
' "${PARAMS_FILE}"; then
  echo "Parameter table validation failed." >&2
  exit 1
fi
parameter_dir="$(cd "$(dirname "${PARAMS_FILE}")" && pwd)"

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
  mkdir -p "${RUN_ROOT}"
  run_root="$(cd "${RUN_ROOT}" && pwd)"
else
  mkdir -p "${SCRATCH}/GPES"
  # Atomically reserve a unique directory even for submissions in the same second.
  run_root="$(mktemp -d "${SCRATCH}/GPES/sweep_$(date +%Y%m%d_%H%M%S)_XXXXXX")"
fi
mkdir -p "${run_root}/slurm"
snapshot="${run_root}/parameter_table.tsv"
if [[ -f "${snapshot}" ]]; then
  if ! cmp -s "${PARAMS_FILE}" "${snapshot}"; then
    echo "RUN_ROOT contains a different parameter table; choose a new RUN_ROOT for a new sweep." >&2
    exit 1
  fi
else
  cp "${PARAMS_FILE}" "${snapshot}"
fi
settings_file="${run_root}/grid_and_seed_settings.txt"
printf -v settings 'NX=%s\nNY=%s\nX_HALF_LENGTH=%s\nY_HALF_LENGTH=%s\nBACKGROUND_WIDTH_X=%s\nBACKGROUND_WIDTH_Y=%s\nSEED_NOISE=%s' \
  "${NX:-1024}" "${NY:-1024}" "${X_HALF_LENGTH:-25}" "${Y_HALF_LENGTH:-25}" \
  "${BACKGROUND_WIDTH_X:-2.5}" "${BACKGROUND_WIDTH_Y:-2.5}" "${SEED_NOISE:-0.001}"
if [[ -f "${settings_file}" ]]; then
  if [[ "$(<"${settings_file}")" != "${settings}" ]]; then
    echo "RUN_ROOT was created with different grid or seed settings; choose a new RUN_ROOT." >&2
    exit 1
  fi
else
  printf '%s\n' "${settings}" > "${settings_file}"
fi

# Peak files are copied at submission time so queued jobs use the same seeds.
mkdir -p "${run_root}/seeds"
row_index=0
while IFS= read -r row || [[ -n "${row}" ]]; do
  [[ -z "${row//[[:space:]]/}" || "${row}" =~ ^[[:space:]]*# ]] && continue
  row_index=$((row_index + 1))
  read -r _a_s _a_dd _particles _dt _duration _wx _wy _wz run_name _rotation seed_file extra <<< "${row}"
  if [[ -n "${extra:-}" ]]; then
    echo "Too many fields in parameter row ${row_index}" >&2
    exit 1
  fi
  if [[ -z "${run_name:-}" ]]; then
    run_name="$(printf 'run_%04d' "${row_index}")"
  fi
  if [[ ! "${run_name}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo "Invalid run name in parameter row ${row_index}: ${run_name}" >&2
    exit 1
  fi
  if [[ -z "${seed_file:-}" || "${seed_file}" == "-" ]]; then
    continue
  fi
  if [[ "${seed_file}" != /* ]]; then
    seed_file="${parameter_dir}/${seed_file}"
  fi
  if [[ ! -f "${seed_file}" ]]; then
    echo "Missing modulation file for ${run_name}: ${seed_file}" >&2
    exit 1
  fi
  submitted_seed="${run_root}/seeds/${run_name}.tsv"
  if [[ -e "${submitted_seed}" ]]; then
    if ! cmp -s "${seed_file}" "${submitted_seed}"; then
      echo "Seed file changed for ${run_name}; choose a new RUN_ROOT for a new seed." >&2
      exit 1
    fi
  else
    cp "${seed_file}" "${submitted_seed}"
  fi
done < "${snapshot}"
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
  "${SCRIPT_DIR}/run.slurm" "${snapshot}")"

echo "Build job: ${build_job}"
echo "Sweep job: ${array_job} (${num_jobs} tasks, at most ${MAX_CONCURRENT} concurrent)"
