cmake_minimum_required(VERSION 3.24)

# Standalone builds consume an installed GPES, with no solver source includes.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(GPESApplications LANGUAGES CXX)
  find_package(gpes 1 CONFIG REQUIRED)
  option(GPES_WARNINGS "Enable compiler warnings" ON)
  option(GPES_BUILD_LEGACY_APP "Build the combined legacy executable" OFF)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
else()
  option(GPES_BUILD_LEGACY_APP "Build the combined legacy executable" ON)
endif()

if(NOT TARGET CLI11::CLI11)
  find_path(GPES_CLI11_INCLUDE_DIR CLI11.hpp
    HINTS "${CMAKE_CURRENT_LIST_DIR}/../include/gpes/third_party/CLI11"
    REQUIRED)
  add_library(CLI11::CLI11 INTERFACE IMPORTED)
  set_target_properties(CLI11::CLI11 PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${GPES_CLI11_INCLUDE_DIR}")
endif()

foreach(program IN ITEMS groundstate bdg)
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} PRIVATE gpes::gpes CLI11::CLI11)
  if(GPES_WARNINGS)
    target_compile_options(${program} PRIVATE
      $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wall;-Wextra;-Wpedantic>)
  endif()
  install(TARGETS ${program} RUNTIME DESTINATION bin)
endforeach()

if(GPES_BUILD_LEGACY_APP)
add_executable(GPESolver
  main.cpp
  cli/cli.cpp
)

# Add executable-only compile flags here.
# Example:
#   -O3
#   -march=native
set(GPES_EXECUTABLE_COMPILE_OPTIONS
    -O3
)

target_compile_options(GPESolver PRIVATE
  ${GPES_EXECUTABLE_COMPILE_OPTIONS}
)

if(GPES_WARNINGS)
  target_compile_options(GPESolver PRIVATE
    $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wall;-Wextra;-Wpedantic>
  )
endif()

target_link_libraries(GPESolver PRIVATE
  CLI11::CLI11
)

if(GPES_ENABLE_PLOTTING)
  target_link_libraries(GPESolver PRIVATE
    gpes::plot
  )
else()
  target_link_libraries(GPESolver PRIVATE
    gpes::gpes
  )
endif()
endif()
