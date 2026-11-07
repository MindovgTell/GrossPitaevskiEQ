add_executable(gpes_numerics_test numerics.cpp)
target_link_libraries(gpes_numerics_test PRIVATE gpes::gpes)
target_compile_options(gpes_numerics_test PRIVATE
  $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wall;-Wextra;-Wpedantic>
)
add_test(NAME gpes_numerics COMMAND gpes_numerics_test)

add_executable(gpes_convergence_test convergence.cpp)
target_link_libraries(gpes_convergence_test PRIVATE gpes::gpes)
add_test(NAME gpes_convergence COMMAND gpes_convergence_test)

add_executable(gpes_rotation_test rotation.cpp)
target_link_libraries(gpes_rotation_test PRIVATE gpes::gpes)
add_test(NAME gpes_rotation COMMAND gpes_rotation_test)

add_executable(gpes_modulations_test modulations.cpp)
target_link_libraries(gpes_modulations_test PRIVATE gpes::gpes)
add_test(NAME gpes_modulations COMMAND gpes_modulations_test)

if(TARGET groundstate AND TARGET bdg)
  add_test(NAME gpes_split_workflow
    COMMAND "${CMAKE_COMMAND}"
      "-DGROUNDSTATE=$<TARGET_FILE:groundstate>"
      "-DBDG=$<TARGET_FILE:bdg>"
      "-DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/split-workflow"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/split_workflow.cmake")
  set_tests_properties(gpes_split_workflow PROPERTIES TIMEOUT 180)
endif()

if(TARGET groundstate)
  add_test(NAME gpes_groundstate_workflow
    COMMAND "${CMAKE_COMMAND}"
      "-DGROUNDSTATE=$<TARGET_FILE:groundstate>"
      "-DREPO_ROOT=${PROJECT_SOURCE_DIR}"
      "-DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/groundstate-workflow"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/groundstate_workflow.cmake")
  set_tests_properties(gpes_groundstate_workflow PROPERTIES TIMEOUT 180)
endif()

add_test(NAME gpes_submit_sweep
  COMMAND "${CMAKE_COMMAND}"
    "-DREPO_ROOT=${PROJECT_SOURCE_DIR}"
    "-DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/submit-sweep"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/submit_sweep.cmake")
set_tests_properties(gpes_submit_sweep PROPERTIES TIMEOUT 30)
