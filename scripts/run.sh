#!/bin/bash 
set -euo pipefail

#Scripts directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#Project directory 
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

#Build directory
GPESOLVER="$PROJECT_ROOT/GPES/build/GPESolver"

echo "Running: $GPESOLVER" 

"$GPESOLVER"
