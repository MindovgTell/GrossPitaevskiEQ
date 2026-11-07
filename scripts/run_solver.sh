#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

EXECUTABLE="${REPO_ROOT}/bin/groundstate"
BUILD_FIRST="0"
RUN_DIR=""
SOLVER_ARGS=()
BUILD_ARGS=()

usage() {
  cat <<EOF
Usage: $0 [options] [-- <groundstate args>]

Runs one local 200x200 ground-state example with a_s=0.009,
a_dd=0.008843, N=80000, and omega_y/omega_x=2.

Options:
  --build                 Build bin/groundstate before running
  --exec PATH             Use a different ground-state executable
  --run-dir PATH          New output directory (default: unique path under runs/)
  --build-arg ARG         Extra argument forwarded to build_local.sh (repeatable)
  -h, --help              Show help

Examples:
  $0 --build
  $0 -- --a-s 0.007 --w-y 3 --max-ground-steps 1000
  $0 -- --resume /path/to/previous/ground_state
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
    --run-dir|--wavefunction-out)
      RUN_DIR="${2:?missing value for --run-dir}"
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

has_solver_option() {
  local option="$1" arg
  if [[ ${#SOLVER_ARGS[@]} -eq 0 ]]; then return 1; fi
  for arg in "${SOLVER_ARGS[@]}"; do
    case "${arg}" in
      "${option}"|"${option}="*) return 0 ;;
    esac
  done
  return 1
}

if has_solver_option --output-dir || has_solver_option --wavefunction-out; then
  echo "Set the output directory with --run-dir, not a forwarded --output-dir." >&2
  exit 1
fi

if [[ "${BUILD_FIRST}" == "1" ]]; then
  if [[ ${#BUILD_ARGS[@]} -gt 0 ]]; then
    "${SCRIPT_DIR}/build_local.sh" "${BUILD_ARGS[@]}"
  else
    "${SCRIPT_DIR}/build_local.sh"
  fi
fi

if [[ ! -x "${EXECUTABLE}" ]]; then
  echo "Executable not found or not executable: ${EXECUTABLE}" >&2
  echo "Run ${SCRIPT_DIR}/build_local.sh first or pass --build." >&2
  exit 1
fi

if [[ -z "${RUN_DIR}" ]]; then
  RUN_DIR="${REPO_ROOT}/runs/local_groundstate_$(date -u +%Y%m%d_%H%M%S)_$$_${RANDOM}"
fi
if [[ -e "${RUN_DIR}" ]]; then
  echo "Output directory already exists: ${RUN_DIR}" >&2
  exit 1
fi
mkdir -p "$(dirname "${RUN_DIR}")"

cmd=("${EXECUTABLE}" --output-dir "${RUN_DIR}")
add_default() {
  if ! has_solver_option "$1"; then
    cmd+=("$1" "$2")
  fi
}

add_default --max-ground-steps 400000
if ! has_solver_option --resume; then
  add_default --nx 200
  add_default --ny 200
  add_default --x-half-length 20
  add_default --y-half-length 20
  add_default --a-s 0.0025
  add_default --a-dd 0.008843
  add_default --num-particles 50000
  add_default --w-x 1
  if ! has_solver_option --lambda; then
    add_default --w-y 5
  fi
  add_default --w-z 8
  add_default --dt 0.0005
fi
if [[ ${#SOLVER_ARGS[@]} -gt 0 ]]; then
  cmd+=("${SOLVER_ARGS[@]}")
fi

echo "Ground-state output: ${RUN_DIR}"
printf 'Command:'
printf ' %q' "${cmd[@]}"
printf '\n'
exec "${cmd[@]}"
