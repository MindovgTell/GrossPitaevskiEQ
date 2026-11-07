cmake_minimum_required(VERSION 3.24)

string(RANDOM LENGTH 10 ALPHABET 0123456789abcdef suffix)
set(root "${TEST_ROOT}/${suffix}")
file(MAKE_DIRECTORY "${root}")
file(WRITE "${root}/peaks.tsv" "# x y sigma_x sigma_y relative_amplitude\n-1 0.2 0.5 0.6 0.8\n0 -0.2 0.5 0.6 0.8\n1 0.2 0.5 0.6 0.8\n")

function(run expected)
  execute_process(COMMAND ${ARGN} RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err TIMEOUT 60)
  if(NOT "${result}" STREQUAL "${expected}")
    message(FATAL_ERROR "Command ${ARGN}\nExpected ${expected}, got ${result}\n${out}\n${err}")
  endif()
endfunction()

function(require_file path)
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "Missing expected output: ${path}")
  endif()
endfunction()

function(require_text path pattern)
  file(READ "${path}" contents)
  if(NOT contents MATCHES "${pattern}")
    message(FATAL_ERROR "Expected ${pattern} in ${path}:\n${contents}")
  endif()
endfunction()

set(args --nx 16 --ny 16 --x-half-length 6 --y-half-length 6
  --a-s 0.001 --a-dd 0 --num-particles 1 --w-x 1 --lambda 0.8 --w-z 5
  --rotation-frequency 0.2 --seed-noise 0 --dt 0.01 --max-ground-steps 1)
run(2 "${GROUNDSTATE}" ${args} --modulation-file "${root}/peaks.tsv" --output-dir "${root}/cli")
require_file("${root}/cli/NOT_CONVERGED")
require_file("${root}/cli/wavefunction/modulations.tsv")
require_file("${root}/cli/images/initial_phase.bmp")
require_file("${root}/cli/images/ground_state_phase.bmp")
require_text("${root}/cli/wavefunction/parameters.txt" "omega_y=1.25")
require_text("${root}/cli/wavefunction/parameters.txt" "rotation_frequency=0.2")
require_text("${root}/cli/ground_state.txt" "modulation_count=3")
require_text("${root}/cli/ground_state.txt" "energy_frame=rotating")
run(2 "${GROUNDSTATE}" --resume "${root}/cli" --max-ground-steps 1 --output-dir "${root}/continued")
require_file("${root}/continued/wavefunction/modulations.tsv")
require_text("${root}/continued/wavefunction/parameters.txt" "rotation_frequency=0.2")

file(MAKE_DIRECTORY "${root}/array/seeds")
file(COPY_FILE "${root}/peaks.tsv" "${root}/array/seeds/testrot.tsv")
file(WRITE "${root}/table.tsv" "0.001 0 1 0.01 0.1 1 1.25 5 testrot 0.2 peaks.tsv\n")
run(2 "${CMAKE_COMMAND}" -E env
  "SLURM_ARRAY_TASK_ID=1" "SLURM_JOB_ID=123" "REPO_ROOT=${REPO_ROOT}"
  "RUN_ROOT=${root}/array" "SCRATCH=${root}"
  "NX=16" "NY=16" "X_HALF_LENGTH=6" "Y_HALF_LENGTH=6"
  "MAX_GROUND_STEPS=1" "SEED_NOISE=0" "RELATIVE_ENERGY_TOL=" "EXECUTABLE=${GROUNDSTATE}"
  bash "${REPO_ROOT}/scripts/run.slurm" "${root}/table.tsv")
require_file("${root}/array/testrot/NOT_CONVERGED")
require_file("${root}/array/testrot/wavefunction/modulations.tsv")
require_text("${root}/array/testrot/wavefunction/parameters.txt" "rotation_frequency=0.2")
file(STRINGS "${root}/array/testrot/ground_state.txt" tolerance_line
  REGEX "^relative_energy_tolerance=")
string(REPLACE "relative_energy_tolerance=" "" tolerance "${tolerance_line}")
if("${tolerance}" LESS 0.0000009 OR "${tolerance}" GREATER 0.0000011)
  message(FATAL_ERROR "Slurm energy tolerance should default to 1e-6; got ${tolerance}")
endif()
message(STATUS "Ground-state CLI and array-task workflow passed: ${root}")
