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

PARAMS_FILE="${1:-${REPO_ROOT}/parameters/sweep_20.tsv}"
EXECUTABLE="${EXECUTABLE:-${REPO_ROOT}/bin/groundstate}"
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
MAX_GROUND_STEPS="${MAX_GROUND_STEPS:-200000}"
MIN_GROUND_DT="${MIN_GROUND_DT:-1e-4}"
CONVERGENCE_CHECK_INTERVAL="${CONVERGENCE_CHECK_INTERVAL:-100}"
CONVERGENCE_CHECKS_REQUIRED="${CONVERGENCE_CHECKS_REQUIRED:-3}"
RELATIVE_ENERGY_TOL="${RELATIVE_ENERGY_TOL:-1e-6}"
ABSOLUTE_ENERGY_TOL="${ABSOLUTE_ENERGY_TOL:-1e-12}"
STATIONARY_RESIDUAL_TOL="${STATIONARY_RESIDUAL_TOL:-1e-4}"
BOUNDARY_DENSITY_TOL="${BOUNDARY_DENSITY_TOL:-1e-8}"
SNAPSHOT_INTERVAL="${SNAPSHOT_INTERVAL:-5000}"
BACKGROUND_WIDTH_X="${BACKGROUND_WIDTH_X:-2.5}"
BACKGROUND_WIDTH_Y="${BACKGROUND_WIDTH_Y:-2.5}"
SEED_NOISE="${SEED_NOISE:-0.001}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

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
  exec "${REPO_ROOT}/scripts/submit_sweep.sh" "${PARAMS_FILE}"
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
  echo "Submit build.slurm first or set EXECUTABLE=/path/to/groundstate." >&2
  exit 1
fi

param_line="$(get_parameter_line "${SLURM_ARRAY_TASK_ID}")"
if [[ -z "${param_line}" ]]; then
  echo "No parameter line for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi

# Parameter-file columns:
# a_s a_dd num_particles dt legacy_duration w_x w_y w_z run_name [Omega] [modulation_file]
# legacy_duration remains in the old sweep table but groundstate uses convergence controls.
read -r a_s a_dd num_particles dt legacy_duration w_x w_y w_z run_name rotation_frequency modulation_file extra <<< "${param_line}"

if [[ -z "${a_s:-}" || -z "${a_dd:-}" || -z "${num_particles:-}" ||
      -z "${dt:-}" || -z "${legacy_duration:-}" ||
      -z "${w_x:-}" || -z "${w_y:-}" || -z "${w_z:-}" ]]; then
  echo "Invalid parameter line ${SLURM_ARRAY_TASK_ID}: ${param_line}" >&2
  echo "Expected: a_s a_dd num_particles dt legacy_duration w_x w_y w_z run_name [Omega] [modulation_file]" >&2
  exit 1
fi

if [[ -n "${extra:-}" ]]; then
  echo "Too many fields on parameter line ${SLURM_ARRAY_TASK_ID}: ${param_line}" >&2
  echo "Expected: a_s a_dd num_particles dt legacy_duration w_x w_y w_z run_name [Omega] [modulation_file]" >&2
  exit 1
fi

if [[ -z "${run_name:-}" ]]; then
  run_name="$(printf "run_%04d" "${SLURM_ARRAY_TASK_ID}")"
fi
if [[ ! "${run_name}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
  echo "Invalid run_name: ${run_name}" >&2
  exit 1
fi
rotation_frequency="${rotation_frequency:-0}"
if [[ -n "${modulation_file:-}" && "${modulation_file}" != "-" ]]; then
  modulation_file="${RUN_ROOT}/seeds/${run_name}.tsv"
  if [[ ! -f "${modulation_file}" ]]; then
    echo "Missing submitted modulation seed: ${modulation_file}" >&2
    exit 1
  fi
else
  modulation_file=""
fi

base_dir="${RUN_ROOT}/${run_name}"
source_dir=""
for candidate in "${base_dir}" "${base_dir}_retry_"*; do
  if [[ -f "${candidate}/COMPLETED" ]]; then
    echo "Ground-state run already complete: ${candidate}"
    exit 0
  fi
  checkpoint="${candidate}/wavefunction/wavefunction.csv"
  if [[ -f "${checkpoint}" && ( -z "${source_dir}" || "${checkpoint}" -nt "${source_dir}/wavefunction/wavefunction.csv" ) ]]; then
    source_dir="${candidate}"
  fi
done

run_dir="${base_dir}"
if [[ -e "${run_dir}" ]]; then
  run_dir="${base_dir}_retry_${SLURM_JOB_ID:-manual}_${SLURM_RESTART_COUNT:-0}"
  attempt=1
  while [[ -e "${run_dir}" ]]; do
    run_dir="${base_dir}_retry_${SLURM_JOB_ID:-manual}_${SLURM_RESTART_COUNT:-0}_${attempt}"
    attempt=$((attempt + 1))
  done
fi

cmd=(
  "${EXECUTABLE}"
  --dt "${dt}"
  --max-ground-steps "${MAX_GROUND_STEPS}"
  --min-ground-dt "${MIN_GROUND_DT}"
  --convergence-check-interval "${CONVERGENCE_CHECK_INTERVAL}"
  --convergence-checks-required "${CONVERGENCE_CHECKS_REQUIRED}"
  --relative-energy-tol "${RELATIVE_ENERGY_TOL}"
  --absolute-energy-tol "${ABSOLUTE_ENERGY_TOL}"
  --stationary-residual-tol "${STATIONARY_RESIDUAL_TOL}"
  --boundary-density-tol "${BOUNDARY_DENSITY_TOL}"
  --snapshot-interval "${SNAPSHOT_INTERVAL}"
  --output-dir "${run_dir}"
)
if [[ -n "${source_dir}" ]]; then
  cmd+=(--resume "${source_dir}")
else
  cmd+=(
    --a-s "${a_s}"
    --a-dd "${a_dd}"
    --num-particles "${num_particles}"
    --w-x "${w_x}"
    --w-y "${w_y}"
    --w-z "${w_z}"
    --rotation-frequency "${rotation_frequency}"
    --nx "${NX}"
    --ny "${NY}"
    --x-half-length "${X_HALF_LENGTH}"
    --y-half-length "${Y_HALF_LENGTH}"
    --background-width-x "${BACKGROUND_WIDTH_X}"
    --background-width-y "${BACKGROUND_WIDTH_Y}"
    --seed-noise "${SEED_NOISE}"
    --random-seed "${SLURM_ARRAY_TASK_ID}"
  )
  if [[ -n "${modulation_file}" ]]; then
    cmd+=(--modulation-file "${modulation_file}")
  fi
fi

echo "SLURM job id: ${SLURM_JOB_ID:-unknown}"
echo "SLURM array task id: ${SLURM_ARRAY_TASK_ID}"
echo "Parameter line: ${param_line}"
echo "Run directory: ${run_dir}"
echo "Resume from: ${source_dir:-fresh seed}"
echo "Grid: ${NX}x${NY} on [-${X_HALF_LENGTH},${X_HALF_LENGTH}) x [-${Y_HALF_LENGTH},${Y_HALF_LENGTH})"
echo "Command: ${cmd[*]}"

exec "${cmd[@]}"
