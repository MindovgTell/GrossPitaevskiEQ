#!/usr/bin/env bash
#SBATCH --job-name=gpes-build
#SBATCH --account=plgquantmol9-cpu
#SBATCH --partition=plgrid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:10:00

set -euo pipefail

# Running this file with bash ignores all #SBATCH directives. Submit it before
# loading the EL9 modules, because the Helios login node still uses RHEL 8.
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "No Slurm allocation detected; submitting build.slurm with sbatch."
  exec sbatch "$0"
fi

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

BUILD_DIR="${BUILD_DIR:-${REPO_ROOT}/build}"
BUILD_TYPE="${BUILD_TYPE:-Release}"
ENABLE_PLOTTING="${ENABLE_PLOTTING:-OFF}"
CLEAN_BUILD="${CLEAN_BUILD:-0}"
ALLOW_FETCHCONTENT="${ALLOW_FETCHCONTENT:-OFF}"
SPECTRA_SOURCE_DIR="${SPECTRA_SOURCE_DIR:-}"
PYTHON_EXECUTABLE="${PYTHON_EXECUTABLE:-}"
PYTHON_VERSION="${PYTHON_VERSION:-}"
GPES_MODULES="${GPES_MODULES:-GCC/14.3.0 CMake/3.31.8 Eigen/3.4.0 FFTW/3.3.10}"

mkdir -p "${REPO_ROOT}/slurm"

cd "${REPO_ROOT}"

if [[ -r /etc/os-release ]]; then
  . /etc/os-release
  OS_MAJOR="${VERSION_ID%%.*}"
  if [[ ! "${OS_MAJOR}" =~ ^[0-9]+$ || "${OS_MAJOR}" -lt 9 ]]; then
    echo "Helios EL9 modules require a Rocky/RHEL 9 compute node, but this job is running on ${PRETTY_NAME:-unknown OS}." >&2
    echo "Submit on Helios with: sbatch scripts/build.slurm" >&2
    exit 1
  fi
fi

echo "Build node: $(hostname)"
echo "Operating system: ${PRETTY_NAME:-unknown}"

if command -v module >/dev/null 2>&1; then
  module purge
  # Helios uses hierarchical Lmod modules. If any module name differs for your
  # selected partition/architecture, check inside the job allocation with:
  #   module spider GCC CMake Eigen FFTW
  module load ${GPES_MODULES}
  module list
else
  echo "Warning: module command not found; using current environment." >&2
fi

if [[ "${CLEAN_BUILD}" == "1" ]]; then
  if [[ "${BUILD_DIR}" != "${REPO_ROOT}/build" ]]; then
    echo "Refusing CLEAN_BUILD for unexpected directory: ${BUILD_DIR}" >&2
    exit 1
  fi
  rm -rf "${BUILD_DIR}"
fi

mkdir -p "${BUILD_DIR}"

CMAKE_ARGS=(
  -S "${REPO_ROOT}"
  -B "${BUILD_DIR}"
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"
  -DGPES_BUILD_APP=ON
  -DGPES_BUILD_TESTS=ON
  -DGPES_ENABLE_PLOTTING="${ENABLE_PLOTTING}"
  -DGPES_ALLOW_FETCHCONTENT="${ALLOW_FETCHCONTENT}"
)

if [[ -n "${SPECTRA_SOURCE_DIR}" ]]; then
  CMAKE_ARGS+=(-DGPES_SPECTRA_SOURCE_DIR="${SPECTRA_SOURCE_DIR}")
fi
if [[ -n "${PYTHON_EXECUTABLE}" ]]; then
  CMAKE_ARGS+=(-DGPES_PYTHON_EXECUTABLE="${PYTHON_EXECUTABLE}")
fi
if [[ -n "${PYTHON_VERSION}" ]]; then
  CMAKE_ARGS+=(-DGPES_PYTHON_VERSION="${PYTHON_VERSION}")
fi

cmake "${CMAKE_ARGS[@]}"

cmake --build "${BUILD_DIR}" --config "${BUILD_TYPE}" --parallel "${SLURM_CPUS_PER_TASK:-4}" \
  --target groundstate gpes_convergence_test gpes_rotation_test gpes_modulations_test
ctest --test-dir "${BUILD_DIR}" --output-on-failure \
  -R '^gpes_(convergence|rotation|modulations|groundstate_workflow|submit_sweep)$'

echo "Ground-state build and numerical tests finished."
echo "Ground-state binary: ${REPO_ROOT}/bin/groundstate"
