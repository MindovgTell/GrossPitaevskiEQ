#!/usr/bin/env bash
set -euo pipefail

# —————————————————————————————————————————
# 1) Load the cluster modules
# —————————————————————————————————————————
# (adjust versions/module‐names as needed — run `module avail` to see what’s installed)
module purge
module load cmake/3.24
module load gcc            # or Intel/oneAPI, etc.
module load python/3.11
module load eigen/3.4
module load fftw

# —————————————————————————————————————————
# 2) Directories
# —————————————————————————————————————————
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
SRC_DIR="$PROJECT_ROOT/GPES"
BUILD_DIR="$SRC_DIR/build"

# make sure build dir exists
mkdir -p "$BUILD_DIR"

# —————————————————————————————————————————
# 3) Configure
# —————————————————————————————————————————
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Configuring GPES in $BUILD_DIR …"
cmake -S "$SRC_DIR" \
      -B "$BUILD_DIR" \
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release} \
      -DCMAKE_CXX_COMPILER="$(which g++)" \
      "$@"

# —————————————————————————————————————————
# 4) Build
# —————————————————————————————————————————
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Building GPES …"
cmake --build "$BUILD_DIR" -- -j"$(nproc)"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] Build complete."