cmake_minimum_required(VERSION 3.24)
project(GPESolver VERSION 1.0.0 LANGUAGES CXX)

# —————————————————————————————————————————
# 1. C++ standard
# —————————————————————————————————————————
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# —————————————————————————————————————————
# 2. Options
# —————————————————————————————————————————
option(GPES_BUILD_APP "Build CLI executable" ON)
option(GPES_ENABLE_PLOTTING "Build plotting module (Python/NumPy + matplotlib-cpp)" ON)
option(GPES_WARNINGS "Enable -Wall -Wextra -Wpedantic" ON)
option(GPES_BUILD_TESTS "Build numerical regression tests" ON)
option(GPES_ALLOW_FETCHCONTENT "Allow CMake to download missing header-only dependencies" OFF)
set(GPES_SPECTRA_SOURCE_DIR "" CACHE PATH "Use an existing Spectra source tree instead of downloading with FetchContent")

# —————————————————————————————————————————
# 3. Platform messages
# —————————————————————————————————————————
if(APPLE)
  message(STATUS "Configuring for macOS")
elseif(UNIX)
  message(STATUS "Configuring for Linux/UNIX")
else()
  message(FATAL_ERROR "Unsupported platform (macOS/Linux only)")
endif()

# —————————————————————————————————————————
# 4. Dependencies (installed packages)
# —————————————————————————————————————————
include(FetchContent)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Eigen
find_package(GPESEigen3 REQUIRED)

# Spectra (header-only; depends on Eigen)
if(GPES_SPECTRA_SOURCE_DIR)
  if(NOT EXISTS "${GPES_SPECTRA_SOURCE_DIR}/include/Spectra/GenEigsSolver.h")
    message(FATAL_ERROR "GPES_SPECTRA_SOURCE_DIR does not look like a Spectra source tree: ${GPES_SPECTRA_SOURCE_DIR}")
  endif()
  set(GPES_SPECTRA_INCLUDE_DIR "${GPES_SPECTRA_SOURCE_DIR}/include")
else()
  find_path(GPES_SPECTRA_INCLUDE_DIR
    NAMES Spectra/GenEigsSolver.h
    HINTS "${CMAKE_CURRENT_SOURCE_DIR}/include/gpes/third_party"
  )
  if(NOT GPES_SPECTRA_INCLUDE_DIR)
    if(NOT GPES_ALLOW_FETCHCONTENT)
      message(FATAL_ERROR "Spectra not found. Load/provide Spectra, set GPES_SPECTRA_SOURCE_DIR, or enable GPES_ALLOW_FETCHCONTENT.")
    endif()
    set(BUILD_TESTS OFF CACHE BOOL "Build dependency tests" FORCE)
    set(BUILD_EXAMPLES OFF CACHE BOOL "Build dependency examples" FORCE)
    FetchContent_Declare(
      spectra
      GIT_REPOSITORY https://github.com/yixuan/spectra.git
      GIT_TAG v1.2.0
      GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(spectra)
    set(GPES_SPECTRA_INCLUDE_DIR "${spectra_SOURCE_DIR}/include")
  endif()
endif()

# FFTW3: search headers/libs directly. Some cluster FFTW CMake package files are
# incomplete, so avoid find_package(FFTW3) unless a target already exists.
find_package(GPESFFTW REQUIRED)
find_package(Threads REQUIRED)


# spdlog: prefer the vendored headers for reproducible offline cluster builds.
set(GPES_BUNDLED_SPDLOG OFF)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/gpes/third_party/spdlog/spdlog.h")
  set(GPES_BUNDLED_SPDLOG ON)
  add_library(spdlog::spdlog_header_only INTERFACE IMPORTED)
  target_compile_definitions(spdlog::spdlog_header_only INTERFACE
    SPDLOG_HEADER_ONLY
    FMT_HEADER_ONLY
  )
  target_include_directories(spdlog::spdlog_header_only INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/include/gpes/third_party"
  )
else()
  find_package(spdlog QUIET)
endif()
if(NOT TARGET spdlog::spdlog_header_only AND NOT TARGET spdlog::spdlog)
  find_path(SPDLOG_INCLUDE_DIR
    NAMES spdlog/spdlog.h
  )
  if(SPDLOG_INCLUDE_DIR)
    add_library(spdlog::spdlog_header_only INTERFACE IMPORTED)
    target_compile_definitions(spdlog::spdlog_header_only INTERFACE SPDLOG_HEADER_ONLY)
    target_include_directories(spdlog::spdlog_header_only INTERFACE "${SPDLOG_INCLUDE_DIR}")
  else()
    if(NOT GPES_ALLOW_FETCHCONTENT)
      message(FATAL_ERROR "spdlog not found. Load/provide spdlog or enable GPES_ALLOW_FETCHCONTENT.")
    endif()
    FetchContent_Declare(
      spdlog
      GIT_REPOSITORY https://github.com/gabime/spdlog.git
      GIT_TAG v1.15.3
      GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(spdlog)
  endif()
endif()

# Only search Python if plotting is enabled
if(GPES_ENABLE_PLOTTING)
  # Python configuration (override with -DGPES_PYTHON_ROOT / -DGPES_PYTHON_EXECUTABLE / -DGPES_PYTHON_VERSION)
  set(GPES_PYTHON_VERSION "" CACHE STRING "Preferred Python version for plotting; empty means use loaded/default Python3")
  set(GPES_PYTHON_ROOT "" CACHE PATH "Preferred Python root directory")
  set(GPES_PYTHON_EXECUTABLE "" CACHE FILEPATH "Preferred Python executable")

  if(APPLE)
    set(Python3_FIND_FRAMEWORK FIRST)
    set(Python3_FIND_STRATEGY LOCATION)
    if(NOT GPES_PYTHON_ROOT AND NOT Python3_ROOT_DIR)
      set(Python3_ROOT_DIR "/Library/Frameworks/Python.framework/Versions/3.11")
    endif()
  endif()

  if(GPES_PYTHON_ROOT)
    set(Python3_ROOT_DIR "${GPES_PYTHON_ROOT}")
  endif()
  if(GPES_PYTHON_EXECUTABLE)
    set(Python3_EXECUTABLE "${GPES_PYTHON_EXECUTABLE}")
  endif()

  if(GPES_PYTHON_VERSION)
    find_package(Python3 ${GPES_PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development NumPy)
  else()
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)
  endif()

  message(STATUS "Python executable: ${Python3_EXECUTABLE}")
  message(STATUS "Python include directory: ${Python3_INCLUDE_DIRS}")
  message(STATUS "Python NumPy include directory: ${Python3_NumPy_INCLUDE_DIRS}")
  message(STATUS "Python libraries: ${Python3_LIBRARIES}")
endif()


# —————————————————————————————————————————
# 5. Header-only CLI11 as an INTERFACE target
# —————————————————————————————————————————
add_library(cli11 INTERFACE)
add_library(CLI11::CLI11 ALIAS cli11)

target_include_directories(cli11 INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/include/gpes/third_party/CLI11
)

# —————————————————————————————————————————
# 6. Core library
# —————————————————————————————————————————
add_library(gpes STATIC
  include/gpes/source/Log/Log.cpp
)

add_library(gpes::gpes ALIAS gpes)

target_include_directories(gpes
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/gpes>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/gpes>"
)
target_compile_features(gpes PUBLIC cxx_std_20)
target_link_libraries(gpes PRIVATE Threads::Threads)

# Link Eigen
target_link_libraries(gpes PUBLIC Eigen3::Eigen)

# Link Spectra
target_include_directories(gpes PUBLIC
    "$<BUILD_INTERFACE:${GPES_SPECTRA_INCLUDE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/gpes/third_party>"
)
# Link FFTW3 (handle both "target" and "variables" styles)
target_link_libraries(gpes PUBLIC FFTW3::fftw3)

# Link spdlog (prefer header_only if available, else compiled target)
if(GPES_BUNDLED_SPDLOG)
  # Give bundled headers an exported target independent of source-tree paths.
  add_library(gpes_logging INTERFACE)
  set_target_properties(gpes_logging PROPERTIES EXPORT_NAME logging)
  target_compile_definitions(gpes_logging INTERFACE SPDLOG_HEADER_ONLY FMT_HEADER_ONLY)
  target_include_directories(gpes_logging INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/gpes/third_party>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/gpes/third_party>"
  )
  target_link_libraries(gpes PUBLIC gpes_logging)
elseif(TARGET spdlog::spdlog_header_only)
  target_link_libraries(gpes PUBLIC spdlog::spdlog_header_only)
  set(GPES_EXTERNAL_SPDLOG_TARGET spdlog::spdlog_header_only)
elseif(TARGET spdlog::spdlog)
  target_link_libraries(gpes PUBLIC spdlog::spdlog)
  set(GPES_EXTERNAL_SPDLOG_TARGET spdlog::spdlog)
else()
  message(FATAL_ERROR "spdlog found but no spdlog::spdlog(_header_only) target available.")
endif()


# --- plotting module (optional) ---
if(GPES_ENABLE_PLOTTING)
  add_library(gpes_plot STATIC
    include/gpes/source/plot/plot.cpp
  )
  add_library(gpes::plot ALIAS gpes_plot)
  set_target_properties(gpes_plot PROPERTIES EXPORT_NAME plot)

  target_include_directories(gpes_plot
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/gpes>"
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/gpes/third_party>"
      "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/gpes>"
      "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/gpes/third_party>"
    PRIVATE
      ${Python3_INCLUDE_DIRS}
      ${Python3_NumPy_INCLUDE_DIRS}
  )

  target_link_libraries(gpes_plot
    PUBLIC
      gpes::gpes
      Python3::Python
      Python3::NumPy
  )

  target_compile_definitions(gpes_plot PUBLIC GPES_HAS_PLOTTING=1)
else()
  # Useful for headers that want to know plotting is unavailable
  target_compile_definitions(gpes PUBLIC GPES_HAS_PLOTTING=0)
endif()

# Warnings
if(GPES_WARNINGS)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(gpes PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endif()

set_target_properties(gpes PROPERTIES POSITION_INDEPENDENT_CODE ON)

include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/InstallGPES.cmake")


# —————————————————————————————————————————
# 7. Output directories
# —————————————————————————————————————————
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
if(CMAKE_CONFIGURATION_TYPES)
  foreach(cfg ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER "${cfg}" cfg_upper)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${cfg_upper} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  endforeach()
endif()

# —————————————————————————————————————————
# 8. CLI app
# —————————————————————————————————————————
if(GPES_BUILD_APP)
  add_subdirectory(app)
endif()

if(GPES_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

# —————————————————————————————————————————
# 8. Final message
# —————————————————————————————————————————
message(STATUS "Configuration complete. Build with:")
