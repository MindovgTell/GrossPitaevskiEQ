cmake_minimum_required(VERSION 3.24)
project(gpes_package_consumer LANGUAGES CXX)

option(GPES_TEST_PLOT "Also check the installed optional plotting component" OFF)
if(GPES_TEST_PLOT)
  find_package(gpes 1 CONFIG REQUIRED COMPONENTS plot)
else()
  find_package(gpes 1 CONFIG REQUIRED)
endif()

add_executable(gpes_consumer main.cpp)
target_link_libraries(gpes_consumer PRIVATE gpes::gpes)
if(GPES_TEST_PLOT)
  target_link_libraries(gpes_consumer PRIVATE gpes::plot)
  target_compile_definitions(gpes_consumer PRIVATE GPES_TEST_PLOT=1)
endif()

enable_testing()
add_test(NAME installed_library_smoke COMMAND gpes_consumer)
if(GPES_TEST_PLOT)
  set_tests_properties(installed_library_smoke PROPERTIES
    ENVIRONMENT "MPLBACKEND=Agg;MPLCONFIGDIR=${CMAKE_CURRENT_BINARY_DIR}/matplotlib")
endif()
