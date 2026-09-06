#!/usr/bin/env bash
#SBATCH --job-name=gpes-run
#SBATCH --account=plgquantmol9-cpu
#SBATCH --partition=plgrid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --requeue

set -euo pipefail

if [[ -z "${REPO_ROOT:-}" ]]; then
  if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO_ROOT="${SLURM_SUBMIT_DIR}"
  else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
  fi
fi

REPO_ROOT="$(cd "${REPO_ROOT}" && pwd)"
if [[ ! -f "${REPO_ROOT}/CMakeLists.txt" && -f "${REPO_ROOT}/../CMakeLists.txt" ]]; then
  REPO_ROOT="$(cd "${REPO_ROOT}/.." && pwd)"
fi
if [[ ! -f "${REPO_ROOT}/CMakeLists.txt" ]]; then
  echo "Could not find repository root. Submit from the repository directory or set REPO_ROOT=/path/to/GrossPitaevski." >&2
  exit 1
fi

PARAMS_FILE="${1:-${REPO_ROOT}/parameters/sweep_40.tsv}"
EXECUTABLE="${EXECUTABLE:-${REPO_ROOT}/bin/GPESolver}"
GPES_MODULES="${GPES_MODULES:-GCC/14.3.0 CMake/3.31.8 Eigen/3.4.0 FFTW/3.3.10}"
if [[ -z "${RUN_ROOT:-}" ]]; then
  if [[ -z "${SCRATCH:-}" ]]; then
    echo "SCRATCH is not set. Set SCRATCH or submit with RUN_ROOT=/path/to/runs." >&2
    exit 1
  fi
  RUN_ROOT="${SCRATCH}/GPES"
fi

NX="${NX:-1024}"
NY="${NY:-1024}"
X_HALF_LENGTH="${X_HALF_LENGTH:-25}"
Y_HALF_LENGTH="${Y_HALF_LENGTH:-25}"
MAX_GROUND_STEPS="${MAX_GROUND_STEPS:-50000}"
MIN_GROUND_DT="${MIN_GROUND_DT:-1e-4}"
CONVERGENCE_CHECK_INTERVAL="${CONVERGENCE_CHECK_INTERVAL:-100}"
SNAPSHOT_INTERVAL="${SNAPSHOT_INTERVAL:-5000}"
BDG_MODES="${BDG_MODES:-24}"
BDG_NCV="${BDG_NCV:-64}"
MAX_CONCURRENT="${MAX_CONCURRENT:-8}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

mkdir -p "${REPO_ROOT}/slurm"

count_parameter_lines() {
  awk '
    /^[[:space:]]*$/ { next }
    /^[[:space:]]*#/ { next }
    { count++ }
    END { print count + 0 }
  ' "$1"
}

get_parameter_line() {
  local line_number="$1"
  awk -v wanted="${line_number}" '
    /^[[:space:]]*$/ { next }
    /^[[:space:]]*#/ { next }
    {
      count++
      if (count == wanted) {
        print
        exit
      }
    }
  ' "${PARAMS_FILE}"
}

if [[ ! -f "${PARAMS_FILE}" ]]; then
  echo "Parameter file not found: ${PARAMS_FILE}" >&2
  exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  num_jobs="$(count_parameter_lines "${PARAMS_FILE}")"
  if [[ "${num_jobs}" -le 0 ]]; then
    echo "No parameter sets found in ${PARAMS_FILE}" >&2
    exit 1
  fi

  echo "Submitting ${num_jobs} GPES jobs from ${PARAMS_FILE}"
  mkdir -p "${RUN_ROOT}/slurm"
  exec sbatch \
    --array="1-${num_jobs}%${MAX_CONCURRENT}" \
    --output="${RUN_ROOT}/slurm/%A_%a.out" \
    --error="${RUN_ROOT}/slurm/%A_%a.err" \
    "$0" "${PARAMS_FILE}"
fi

if [[ -r /etc/os-release ]]; then
  . /etc/os-release
  OS_MAJOR="${VERSION_ID%%.*}"
  if [[ ! "${OS_MAJOR}" =~ ^[0-9]+$ || "${OS_MAJOR}" -lt 9 ]]; then
    echo "Helios EL9 modules require a Rocky/RHEL 9 compute node, but this job is running on ${PRETTY_NAME:-unknown OS}." >&2
    echo "Submit on Helios with: sbatch scripts/run.slurm" >&2
    exit 1
  fi
fi

if command -v module >/dev/null 2>&1; then
  module purge
  # Keep runtime library paths consistent with build.slurm.
  module load ${GPES_MODULES}
  module list
else
  echo "Warning: module command not found; using current environment." >&2
fi

if [[ ! -x "${EXECUTABLE}" ]]; then
  echo "Executable not found or not executable: ${EXECUTABLE}" >&2
  echo "Submit build.slurm first or set EXECUTABLE=/path/to/GPESolver." >&2
  exit 1
fi

param_line="$(get_parameter_line "${SLURM_ARRAY_TASK_ID}")"
if [[ -z "${param_line}" ]]; then
  echo "No parameter line for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi

# Parameter-file columns:
# a_s a_dd num_particles dt duration w_x w_y w_z [run_name]
read -r a_s a_dd num_particles dt duration w_x w_y w_z run_name extra <<< "${param_line}"

if [[ -z "${a_s:-}" || -z "${a_dd:-}" || -z "${num_particles:-}" ||
      -z "${dt:-}" || -z "${duration:-}" ||
      -z "${w_x:-}" || -z "${w_y:-}" || -z "${w_z:-}" ]]; then
  echo "Invalid parameter line ${SLURM_ARRAY_TASK_ID}: ${param_line}" >&2
  echo "Expected: a_s a_dd num_particles dt duration w_x w_y w_z [run_name]" >&2
  exit 1
fi

if [[ -n "${extra:-}" ]]; then
  echo "Too many fields on parameter line ${SLURM_ARRAY_TASK_ID}: ${param_line}" >&2
  echo "Expected: a_s a_dd num_particles dt duration w_x w_y w_z [run_name]" >&2
  exit 1
fi

if [[ -z "${run_name:-}" ]]; then
  run_name="$(printf "run_%04d" "${SLURM_ARRAY_TASK_ID}")"
fi

run_dir="${RUN_ROOT}/${run_name}"
log_dir="${run_dir}/log"

mkdir -p "${run_dir}" "${log_dir}"

if [[ -f "${run_dir}/COMPLETED" ]]; then
  echo "Run already complete: ${run_dir}"
  exit 0
fi

cmd=(
  "${EXECUTABLE}"
  --cli-only
  --a-s "${a_s}"
  --a-dd "${a_dd}"
  --num-particles "${num_particles}"
  --dt "${dt}"
  --duration "${duration}"
  --w-x "${w_x}"
  --w-y "${w_y}"
  --w-z "${w_z}"
  --nx "${NX}"
  --ny "${NY}"
  --x-half-length "${X_HALF_LENGTH}"
  --y-half-length "${Y_HALF_LENGTH}"
  --max-ground-steps "${MAX_GROUND_STEPS}"
  --min-ground-dt "${MIN_GROUND_DT}"
  --convergence-check-interval "${CONVERGENCE_CHECK_INTERVAL}"
  --snapshot-interval "${SNAPSHOT_INTERVAL}"
  --random-seed "${SLURM_ARRAY_TASK_ID}"
  --bdg-modes "${BDG_MODES}"
  --bdg-ncv "${BDG_NCV}"
  --log-dir "${log_dir}"
  --wavefunction-out "${run_dir}"
)

echo "SLURM job id: ${SLURM_JOB_ID:-unknown}"
echo "SLURM array task id: ${SLURM_ARRAY_TASK_ID}"
echo "Parameter line: ${param_line}"
echo "Run directory: ${run_dir}"
echo "Grid: ${NX}x${NY} on [-${X_HALF_LENGTH},${X_HALF_LENGTH}) x [-${Y_HALF_LENGTH},${Y_HALF_LENGTH})"
echo "Command: ${cmd[*]}"

exec "${cmd[@]}"
