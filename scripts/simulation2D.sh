#!/bin/bash 
set -euo pipefail

#Scripts directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#Project directory 
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

#Build directory
GPESOLVER="$PROJECT_ROOT/GPES/build/GPESolver"

OUT_DIR="$PROJECT_ROOT/res/2D/Fin"


# set the parameters for simulations
GridSize=(100 200 300 400 500 600 700 800 900 1000)
GridSizeX='300'
GridSizeY='300'
TimeStep='0.001'
StartX='-30'
StartY='-30'
NumberOfMolecules='1000'
e_dd_list=(0.2)
a_s_list=(0.1)
# e_dd='0.9'
# a_s='0.1'
sigma_x='5'
sigma_y='5'
omega_x='1'
omega_y='2'

mkdir -p "$OUT_DIR"

for a_s in "${a_s_list[@]}"; do
  for e_dd in "${e_dd_list[@]}"; do

    # construct output filename embedding the params
    RUN_DIR="$OUT_DIR/finstate_a_s${a_s}_e_dd${e_dd}_omega_y2"
    mkdir -p "$RUN_DIR"

    echo "Running: $GPESOLVER \
      $GridSizeX $GridSizeY $TimeStep \
      $StartX $StartY $NumberOfMolecules \
      $e_dd $a_s $sigma_x $sigma_y $omega_x $omega_y $RUN_DIR"

    # invoke solver (sequentially), redirecting its CSV to OUT_FILE
    "$GPESOLVER" \
      "$GridSizeX" "$GridSizeY" "$TimeStep" \
      "$StartX"  "$StartY"  "$NumberOfMolecules" \
      "$e_dd"    "$a_s"     "$sigma_x" \
      "$sigma_y" "$omega_x" "$omega_y" "$RUN_DIR"

  done
done