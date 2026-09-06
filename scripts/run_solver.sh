#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

EXECUTABLE="${REPO_ROOT}/bin/GPESolver"
BUILD_FIRST="0"
RUN_DIR=""
LOG_DIR=""
WAVEFUNCTION_OUT=""
SOLVER_ARGS=()
BUILD_ARGS=()

usage() {
  cat <<EOF
Usage: $0 [options] [-- <solver args>]

Options:
  --build                 Run scripts/build_local.sh before starting the solver
  --exec PATH             Path to solver executable
  --run-dir PATH          Directory for logs and outputs
  --log-dir PATH          Explicit log directory passed to the solver
  --wavefunction-out PATH Explicit run output directory passed to the solver
  --build-arg ARG         Extra argument forwarded to build_local.sh (repeatable)
  -h, --help              Show help

Examples:
  $0 --build -- --a-s 0.00607 --a-dd 0.00882
  $0 --run-dir ./runs/test -- --dt 0.0005 --duration 0.2
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --build)
      BUILD_FIRST="1"
      shift
      ;;
    --exec)
      EXECUTABLE="${2:?missing value for --exec}"
      shift 2
      ;;
    --run-dir)
      RUN_DIR="${2:?missing value for --run-dir}"
      shift 2
      ;;
    --log-dir)
      LOG_DIR="${2:?missing value for --log-dir}"
      shift 2
      ;;
    --wavefunction-out)
      WAVEFUNCTION_OUT="${2:?missing value for --wavefunction-out}"
      shift 2
      ;;
    --build-arg)
      BUILD_ARGS+=("${2:?missing value for --build-arg}")
      shift 2
      ;;
    --)
      shift
      SOLVER_ARGS+=("$@")
      break
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      SOLVER_ARGS+=("$1")
      shift
      ;;
  esac
done

if [[ "${BUILD_FIRST}" == "1" ]]; then
  "${SCRIPT_DIR}/build_local.sh" "${BUILD_ARGS[@]}"
fi

if [[ -z "${RUN_DIR}" ]]; then
  timestamp="$(date +%Y%m%d_%H%M%S)"
  RUN_DIR="${REPO_ROOT}/runs/${timestamp}"
fi

if [[ -z "${LOG_DIR}" ]]; then
  LOG_DIR="${RUN_DIR}/logs"
fi

if [[ -z "${WAVEFUNCTION_OUT}" ]]; then
  WAVEFUNCTION_OUT="${RUN_DIR}"
fi

mkdir -p "${LOG_DIR}"
mkdir -p "$(dirname "${WAVEFUNCTION_OUT}")"

if [[ ! -x "${EXECUTABLE}" ]]; then
  echo "Executable not found or not executable: ${EXECUTABLE}" >&2
  echo "Run ${SCRIPT_DIR}/build_local.sh first or pass --build." >&2
  exit 1
fi

cmd=("${EXECUTABLE}" "--log-dir" "${LOG_DIR}" "--wavefunction-out" "${WAVEFUNCTION_OUT}")

if [[ ${#SOLVER_ARGS[@]} -gt 0 ]]; then
  cmd+=("${SOLVER_ARGS[@]}")
fi

echo "Run directory: ${RUN_DIR}"
echo "Log directory: ${LOG_DIR}"
echo "Wavefunction:  ${WAVEFUNCTION_OUT}"
echo "Command: ${cmd[*]}"

exec "${cmd[@]}"
