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

# FFTW3
# Some distros provide FFTW3::fftw3, others provide variables.
find_package(FFTW3 REQUIRED)

# spdlog (installed)
# On most systems this provides spdlog::spdlog or spdlog::spdlog_header_only
find_package(spdlog REQUIRED)

# Only search Python if plotting is enabled
if(GPES_ENABLE_PLOTTING)
  find_package(Python3 3.11 REQUIRED COMPONENTS Interpreter Development NumPy)
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
    ${CMAKE_CURRENT_SOURCE_DIR}/include
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
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/third_party/matplotlib-cpp
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
# 7. CLI app
# —————————————————————————————————————————
if(GPES_BUILD_APP)
  add_subdirectory(app)
endif()

# —————————————————————————————————————————
# 8. Final message
# —————————————————————————————————————————
message(STATUS "Configuration complete. Build with:")


# cmake_minimum_required(VERSION 3.24)
# project(GPESolver VERSION 1.0.0 LANGUAGES CXX)

# # —————————————————————————————————————————
# # 1. C++ standard
# # —————————————————————————————————————————
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

# # —————————————————————————————————————————
# # 2. Platform-specific configuration
# # —————————————————————————————————————————
# if(APPLE)
#   message(STATUS "Configuring for macOS (Homebrew)")
  
#   # Python3 (Homebrew)
#   set(Python3_ROOT_DIR "/Library/Frameworks/Python.framework/Versions/3.11")
#   set(Python3_EXECUTABLE "${Python3_ROOT_DIR}/bin/python3.11")
#   find_package(Python3 3.11 REQUIRED COMPONENTS Interpreter Development NumPy)
  
#   # Eigen3 (Homebrew)
#   find_path(EIGEN3_INCLUDE_DIR NAMES Eigen/Core
#       PATHS /opt/homebrew/include /usr/local/include /usr/include
#       PATH_SUFFIXES eigen3)
#   if(NOT EIGEN3_INCLUDE_DIR)
#       message(FATAL_ERROR "Eigen3 not found on macOS. Please brew install eigen")
#   endif()
  
#   # FFTW3 (Homebrew)
#   find_path(FFTW3_INCLUDE_DIR fftw3.h
#       PATHS /opt/homebrew/include /usr/local/include /usr/include)
#   find_library(FFTW3_LIBRARY
#       NAMES fftw3
#       PATHS /opt/homebrew/lib /usr/local/lib /usr/lib)
#   if(NOT FFTW3_INCLUDE_DIR OR NOT FFTW3_LIBRARY)
#       message(FATAL_ERROR "FFTW3 not found on macOS. Please brew install fftw")
#   endif()
#   set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
#   set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})
  
# else()
#   message(STATUS "Configuring for PL-Grid (Linux)")
  
#   # Python3
#   find_package(Python3 3.11 REQUIRED COMPONENTS Interpreter Development NumPy)
  
#   # Eigen3 (header-only)
#   find_package(Eigen3 3.3 REQUIRED NO_MODULE)
  
#   # FFTW3
#   find_package(FFTW3 REQUIRED)
  
#   # Threads (for fftw3_threads if needed)
#   find_package(Threads)
# endif()

# # —————————————————————————————————————————
# # 3. Executable and linking
# # —————————————————————————————————————————
# add_executable(${PROJECT_NAME} ../app/main.cpp)

# target_include_directories(${PROJECT_NAME} PRIVATE
#   ${Python3_INCLUDE_DIRS}
#   ${Python3_NumPy_INCLUDE_DIRS}
#   ${CMAKE_CURRENT_SOURCE_DIR}            # for matplotlibcpp.h
#   $<$<BOOL:${EIGEN3_INCLUDE_DIR}>:${EIGEN3_INCLUDE_DIR}>  # macOS Eigen or PL-Grid Eigen3::Eigen
#   $<$<TARGET_EXISTS:Eigen3::Eigen>:${Eigen3_INCLUDE_DIR}> # PL-Grid Eigen
#   $<$<BOOL:${FFTW3_INCLUDE_DIRS}>:${FFTW3_INCLUDE_DIRS}>  # macOS FFTW
#   $<$<TARGET_EXISTS:FFTW3::fftw3>:${FFTW3_INCLUDE_DIRS}>  # PL-Grid FFTW
# )

# target_link_libraries(${PROJECT_NAME} PRIVATE
#   Python3::Python
#   Python3::NumPy
#   $<$<TARGET_EXISTS:Eigen3::Eigen>:Eigen3::Eigen>
#   $<$<TARGET_EXISTS:FFTW3::fftw3>:FFTW3::fftw3>
#   $<$<BOOL:${FFTW3_LIBRARIES}>:${FFTW3_LIBRARIES}>  # macOS FFTW lib
# #   Threads::Threads
# )

# # —————————————————————————————————————————
# # 4. Final messages
# # —————————————————————————————————————————
# message(STATUS "Configuration complete. Build with: mkdir build && cd build && cmake .. && make -j")