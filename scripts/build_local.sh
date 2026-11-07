#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BUILD_DIR="${REPO_ROOT}/build"
BUILD_TYPE="Release"
ENABLE_PLOTTING="OFF"
CLEAN_BUILD="0"
EXTRA_CMAKE_ARGS=()

usage() {
  cat <<EOF
Usage: $0 [options] [-- <extra cmake configure args>]

Options:
  --build-dir PATH        Build directory (default: ${BUILD_DIR})
  --build-type TYPE       CMake build type (default: ${BUILD_TYPE})
  --plotting              Configure with -DGPES_ENABLE_PLOTTING=ON (default: OFF)
  --no-plotting           Keep plotting disabled
  --clean                 Remove the build directory before configuring
  -h, --help              Show help

Examples:
  $0
  $0 --build-type Debug
  $0 --plotting -- -DGPES_PYTHON_EXECUTABLE=/path/to/python3
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --build-dir)
      BUILD_DIR="${2:?missing value for --build-dir}"
      shift 2
      ;;
    --build-type)
      BUILD_TYPE="${2:?missing value for --build-type}"
      shift 2
      ;;
    --plotting)
      ENABLE_PLOTTING="ON"
      shift
      ;;
    --no-plotting)
      ENABLE_PLOTTING="OFF"
      shift
      ;;
    --clean)
      CLEAN_BUILD="1"
      shift
      ;;
    --)
      shift
      EXTRA_CMAKE_ARGS+=("$@")
      break
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ "${CLEAN_BUILD}" == "1" ]]; then
  if [[ "${BUILD_DIR}" != "${REPO_ROOT}/build" ]]; then
    echo "Refusing --clean for unexpected directory: ${BUILD_DIR}" >&2
    exit 1
  fi
  rm -rf "${BUILD_DIR}"
fi

mkdir -p "${BUILD_DIR}"

cmake_cmd=(
  cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}"
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"
  -DGPES_BUILD_APP=ON
  -DGPES_ENABLE_PLOTTING="${ENABLE_PLOTTING}"
)

if [[ ${#EXTRA_CMAKE_ARGS[@]} -gt 0 ]]; then
  cmake_cmd+=("${EXTRA_CMAKE_ARGS[@]}")
fi

"${cmake_cmd[@]}"

cmake --build "${BUILD_DIR}" --config "${BUILD_TYPE}" --target groundstate

echo "Ground-state build finished."
echo "Binary: ${REPO_ROOT}/bin/groundstate"
