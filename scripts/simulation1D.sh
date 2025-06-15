#!/bin/bash 
set -euo pipefail

#Scripts directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#Project directory 
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

#Build directory
GPESOLVER="$PROJECT_ROOT/build/GPESolver"

OUT_DIR="$PROJECT_ROOT/res/2D/Fin"

# set the parameters for simulations
GridSizeX='1000'
TimeStep='0.0001'
StartX='-30'
NumberOfMolecules='1000'
e_dd='0.1'
a_s='0.1'
sigma_x='5'
omega_x='1'

OUT_FILE="$OUT_DIR/finstate_a_s${a_s}_e_dd${e_dd}.csv"
mkdir -p "$OUT_DIR"

echo "Running: $GPESOLVER \
  $GridSizeX $TimeStep \
  $StartX $NumberOfMolecules \
  $e_dd $a_s $sigma_x $omega_x"

$GPESOLVER \
  "$GridSizeX" "$TimeStep" \
  "$StartX" "$NumberOfMolecules" \
  "$e_dd"    "$a_s"     "$sigma_x" \
  "$omega_x"