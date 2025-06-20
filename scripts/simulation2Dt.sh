#!/bin/bash 
set -euo pipefail

#Scripts directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#Project directory 
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

#Build directory
GPESOLVER="$PROJECT_ROOT/GPES/build/GPESolver"

OUT_DIR="$PROJECT_ROOT/res/2D/FinTest"


# set the parameters for simulations
GridSizeX='100'
GridSizeY='100'
TimeStep='0.001'
StartX='-30'
StartY='-30'
NumberOfMolecules='1000'
# e_dd_list=(1.0 1.15 1.25 1.35 1.4)
# a_s_list=(0.1)
e_dd='0.9'
a_s='0.1'
sigma_x='5'
sigma_y='5'
omega_x='1'
omega_y='1'

OUT_FILE="$OUT_DIR/finstate_a_s${a_s}_e_dd${e_dd}_5"


mkdir -p "$OUT_DIR"

echo "Running: $GPESOLVER \
  $GridSizeX $GridSizeY $TimeStep \
  $StartX $StartY $NumberOfMolecules \
  $e_dd $a_s $sigma_x $sigma_y $omega_x $omega_y $OUT_FILE"

$GPESOLVER \
  "$GridSizeX" "$GridSizeY" "$TimeStep" \
  "$StartX"  "$StartY"  "$NumberOfMolecules" \
  "$e_dd"    "$a_s"     "$sigma_x" \
  "$sigma_y" "$omega_x" "$omega_y" "$OUT_FILE"
