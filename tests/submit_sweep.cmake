cmake_minimum_required(VERSION 3.24)

string(RANDOM LENGTH 10 ALPHABET 0123456789abcdef suffix)
set(root "${TEST_ROOT}/${suffix}")
file(MAKE_DIRECTORY "${root}/mock")
file(WRITE "${root}/mock/sbatch" [=[#!/usr/bin/env bash
printf '%s\n' "$*" >> "$SBATCH_LOG"
if [[ "$*" == *build.slurm* ]]; then
  printf '101\n'
else
  printf '102\n'
fi
]=])
file(CHMOD "${root}/mock/sbatch" PERMISSIONS
  OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(WRITE "${root}/seed.tsv" "0 0 0.5 0.5 0.8\n")
file(WRITE "${root}/table.tsv" "0.001 0.008843 10 0.01 0.1 1 2 8 trial 0.2 seed.tsv\n")

set(ENV{RUN_ROOT} "")
set(ENV{MAX_CONCURRENT} "")

function(run expected)
  execute_process(COMMAND ${ARGN} RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err TIMEOUT 20)
  if(NOT "${result}" STREQUAL "${expected}")
    message(FATAL_ERROR "Command ${ARGN}\nExpected ${expected}, got ${result}\n${out}\n${err}")
  endif()
  set(last_output "${out}" PARENT_SCOPE)
endfunction()

set(command "${CMAKE_COMMAND}" -E env
  "SCRATCH=${root}" "PATH=${root}/mock:$ENV{PATH}" "SBATCH_LOG=${root}/sbatch.log"
  bash "${REPO_ROOT}/scripts/submit_sweep.sh" "${root}/table.tsv")
run(0 ${command})
if(NOT last_output MATCHES "Results:   ([^\n]+)")
  message(FATAL_ERROR "Submission did not print its result directory: ${last_output}")
endif()
set(first "${CMAKE_MATCH_1}")
foreach(path "${first}/parameter_table.tsv" "${first}/grid_and_seed_settings.txt"
             "${first}/seeds/trial.tsv" "${first}/slurm")
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "Submission did not create ${path}")
  endif()
endforeach()
file(READ "${first}/grid_and_seed_settings.txt" settings)
if(NOT settings MATCHES "NX=1024" OR NOT settings MATCHES "NY=1024")
  message(FATAL_ERROR "Wrong grid settings snapshot: ${settings}")
endif()
file(READ "${root}/sbatch.log" calls)
if(NOT calls MATCHES "--dependency=afterok:101" OR NOT calls MATCHES "--array=1-1%8")
  message(FATAL_ERROR "Wrong Slurm dependency or array: ${calls}")
endif()

run(0 ${command})
if(NOT last_output MATCHES "Results:   ([^\n]+)")
  message(FATAL_ERROR "Second submission did not print its result directory: ${last_output}")
endif()
set(second "${CMAKE_MATCH_1}")
if(first STREQUAL second)
  message(FATAL_ERROR "Two default submissions reused ${first}")
endif()

run(0 "${CMAKE_COMMAND}" -E env "RUN_ROOT=${first}" ${command})
run(1 "${CMAKE_COMMAND}" -E env "RUN_ROOT=${first}" "NX=32" ${command})
file(WRITE "${root}/seed.tsv" "0 0 0.5 0.5 0.9\n")
run(1 "${CMAKE_COMMAND}" -E env "RUN_ROOT=${first}" ${command})

file(STRINGS "${REPO_ROOT}/parameters/sweep_80.tsv" default_rows REGEX "^[0-9]")
list(LENGTH default_rows default_count)
if(NOT default_count EQUAL 80)
  message(FATAL_ERROR "Default sweep has ${default_count} rows, expected 80")
endif()
set(scattering_lengths)
set(ratios)
foreach(row IN LISTS default_rows)
  string(REGEX REPLACE "[ \t]+" ";" fields "${row}")
  list(GET fields 0 a_s)
  list(GET fields 1 a_dd)
  list(GET fields 6 omega_y)
  list(GET fields 7 omega_z)
  if(NOT a_dd STREQUAL "0.008843" OR NOT omega_y MATCHES "^[1-5]$" OR NOT omega_z STREQUAL "8")
    message(FATAL_ERROR "Invalid default sweep row: ${row}")
  endif()
  list(APPEND scattering_lengths "${a_s}")
  list(APPEND ratios "${omega_y}")
endforeach()
list(REMOVE_DUPLICATES scattering_lengths)
list(REMOVE_DUPLICATES ratios)
list(LENGTH scattering_lengths a_count)
list(LENGTH ratios ratio_count)
if(NOT a_count EQUAL 16 OR NOT ratio_count EQUAL 5)
  message(FATAL_ERROR "Expected 16 scattering lengths x 5 ratios; got ${a_count} x ${ratio_count}")
endif()
run(0 "${CMAKE_COMMAND}" -E env
  "SCRATCH=${root}" "PATH=${root}/mock:$ENV{PATH}" "SBATCH_LOG=${root}/default_sbatch.log"
  bash "${REPO_ROOT}/scripts/submit_sweep.sh")
file(READ "${root}/default_sbatch.log" default_calls)
if(NOT default_calls MATCHES "--array=1-80%8")
  message(FATAL_ERROR "Default submission did not schedule 80 tasks: ${default_calls}")
endif()
message(STATUS "Submission snapshot and retry guards passed: ${root}")
