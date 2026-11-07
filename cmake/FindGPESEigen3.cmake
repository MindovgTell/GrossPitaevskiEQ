if(TARGET Eigen3::Eigen)
  set(GPESEigen3_FOUND TRUE)
  return()
endif()

find_package(Eigen3 3.3 QUIET CONFIG)
if(TARGET Eigen3::Eigen)
  set(GPESEigen3_FOUND TRUE)
  return()
endif()

find_path(EIGEN3_INCLUDE_DIR NAMES Eigen/Core PATH_SUFFIXES eigen3 Eigen3)
set(GPESEigen3_VERSION "")
if(EIGEN3_INCLUDE_DIR)
  file(STRINGS "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h"
    _gpes_eigen_versions REGEX "^#define EIGEN_(WORLD|MAJOR|MINOR)_VERSION [0-9]+")
  foreach(_part WORLD MAJOR MINOR)
    string(REGEX MATCH "EIGEN_${_part}_VERSION ([0-9]+)" _match "${_gpes_eigen_versions}")
    set(_gpes_eigen_${_part} "${CMAKE_MATCH_1}")
  endforeach()
  set(GPESEigen3_VERSION "${_gpes_eigen_WORLD}.${_gpes_eigen_MAJOR}.${_gpes_eigen_MINOR}")
endif()

set(_gpes_eigen_version_ok FALSE)
if(GPESEigen3_VERSION VERSION_GREATER_EQUAL 3.3)
  set(_gpes_eigen_version_ok TRUE)
endif()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GPESEigen3
  REQUIRED_VARS EIGEN3_INCLUDE_DIR _gpes_eigen_version_ok
  VERSION_VAR GPESEigen3_VERSION
  REASON_FAILURE_MESSAGE "Provide Eigen >= 3.3 through a module, CMAKE_PREFIX_PATH, or EIGEN3_INCLUDE_DIR.")
if(GPESEigen3_FOUND)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
endif()
mark_as_advanced(EIGEN3_INCLUDE_DIR)
