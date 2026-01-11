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

# Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)


if(APPLE)
  message(STATUS "Configuring for macOS (Homebrew)")
  
  # FFTW3 (Homebrew)
  find_path(FFTW3_INCLUDE_DIR fftw3.h
      PATHS /opt/homebrew/include /usr/local/include /usr/include)
  find_library(FFTW3_LIBRARY
      NAMES fftw3
      PATHS /opt/homebrew/lib /usr/local/lib /usr/lib)
  if(NOT FFTW3_INCLUDE_DIR OR NOT FFTW3_LIBRARY)
      message(FATAL_ERROR "FFTW3 not found on macOS. Please brew install fftw")
  endif()
  set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
  set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})
else()
  message(STATUS "Configuring for Linux/UNIX")
  
  # FFTW3 (Linux)
  find_package(FFTW3 REQUIRED)
endif()


# spdlog (installed)
# On most systems this provides spdlog::spdlog or spdlog::spdlog_header_only
find_package(spdlog REQUIRED)

# Only search Python if plotting is enabled
if(GPES_ENABLE_PLOTTING)
  # Python configuration (override with -DGPES_PYTHON_ROOT / -DGPES_PYTHON_EXECUTABLE / -DGPES_PYTHON_VERSION)
  set(GPES_PYTHON_VERSION "3.11" CACHE STRING "Preferred Python version for plotting")
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

  find_package(Python3 ${GPES_PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development NumPy)

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
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/CLI11
)

# —————————————————————————————————————————
# 6. Core library
# —————————————————————————————————————————
add_library(gpes STATIC
  source/Log/log.cpp
)

add_library(gpes::gpes ALIAS gpes)

target_include_directories(gpes
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include/gpes
)

# Link Eigen
target_link_libraries(gpes PUBLIC Eigen3::Eigen)

# Link FFTW3 (handle both "target" and "variables" styles)
if(TARGET FFTW3::fftw3)
  target_link_libraries(gpes PUBLIC FFTW3::fftw3)
else()
  # Common variables exported by FindFFTW3.cmake:
  # FFTW3_LIBRARIES, FFTW3_INCLUDE_DIRS
  if(FFTW3_INCLUDE_DIRS)
    target_include_directories(gpes PUBLIC ${FFTW3_INCLUDE_DIRS})
  endif()
  if(FFTW3_LIBRARIES)
    target_link_libraries(gpes PUBLIC ${FFTW3_LIBRARIES})
  else()
    message(FATAL_ERROR "FFTW3 found but no FFTW3::fftw3 target and no FFTW3_LIBRARIES provided.")
  endif()
endif()

# Link spdlog (prefer header_only if available, else compiled target)
if(TARGET spdlog::spdlog_header_only)
  target_link_libraries(gpes PUBLIC spdlog::spdlog_header_only)
elseif(TARGET spdlog::spdlog)
  target_link_libraries(gpes PUBLIC spdlog::spdlog)
else()
  message(FATAL_ERROR "spdlog found but no spdlog::spdlog(_header_only) target available.")
endif()


# --- plotting module (optional) ---
if(GPES_ENABLE_PLOTTING)
  add_library(gpes_plot STATIC
    source/plot/plot.cpp
  )
  add_library(gpes::plot ALIAS gpes_plot)

  target_include_directories(gpes_plot
    PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR}/include
      ${CMAKE_CURRENT_SOURCE_DIR}/third_party
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

# —————————————————————————————————————————
# 8. Final message
# —————————————————————————————————————————
message(STATUS "Configuration complete. Build with:")
