#!/bin/bash 
set -euo pipefail


#Scripts directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#Project directory 
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

#Source files directory
SRC_DIR="$PROJECT_ROOT/GPES"
#Build directory
# BUILD_DIR="$PROJECT_ROOT/GPES/build"
BUILD_DIR="$PROJECT_ROOT/build"

mkdir -p "$BUILD_DIR"
#Configure & build
echo "Configuring GPES in $BUILD_DIR …"
cmake -S "$PROJECT_ROOT" -B "$BUILD_DIR" "$@"

echo "Building GPES …"
cmake --build "$BUILD_DIR" -- -j"$(nproc)"

echo "Configureation is Done!"

