# cmake_minimum_required(VERSION 3.27.5)

# set(PROJECT_NAME GPESolver)
# project(${PROJECT_NAME} VERSION 1.0.0)

# # Set C++ standard
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED True)

# # Python configuration
# set(Python3_ROOT_DIR "/Library/Frameworks/Python.framework/Versions/3.11")
# set(Python3_EXECUTABLE "/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11")

# find_package(Python3 3.11 REQUIRED COMPONENTS Interpreter Development NumPy)

# message(STATUS "Python executable: ${Python3_EXECUTABLE}")
# message(STATUS "Python include directory: ${Python3_INCLUDE_DIRS}")
# message(STATUS "Python NumPy include directory: ${Python3_NumPy_INCLUDE_DIRS}")
# message(STATUS "Python libraries: ${Python3_LIBRARIES}")

# # Eigen configuration
# find_path(EIGEN3_INCLUDE_DIR NAMES Eigen/Core
#     PATHS /opt/homebrew/include /usr/include /usr/local/include
#     PATH_SUFFIXES eigen3)

# if(NOT EIGEN3_INCLUDE_DIR)
#     message(FATAL_ERROR "Can't find Eigen")
# endif()

# message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}")

# find_path(FFTW3_INCLUDE_DIR fftw3.h
#     PATHS /opt/homebrew/include /usr/local/include /usr/include
# )

# find_library(FFTW3_LIBRARY
#     NAMES fftw3
#     PATHS /opt/homebrew/lib /usr/local/lib /usr/lib
# )

# if(FFTW3_INCLUDE_DIR AND FFTW3_LIBRARY)
#     message(STATUS "FFTW3 found manually")
#     set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
#     set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})
# else()
#     message(FATAL_ERROR "FFTW3 not found - please install via 'brew install fftw'")
# endif()
# # endif()

# # Create executable
# add_executable(${PROJECT_NAME} main.cpp)

# # Include directories
# target_include_directories(${PROJECT_NAME} PRIVATE
#     ${Python3_INCLUDE_DIRS}        # Python 3.11
#     ${Python3_NumPy_INCLUDE_DIRS}  # NumPy
#     ${EIGEN3_INCLUDE_DIR}          # Eigen
#     ${CMAKE_SOURCE_DIR}            # matplotlibcpp.h
#     ${FFTW3_INCLUDE_DIRS}          # FFTW3
# )

# # Link libraries
# target_link_libraries(${PROJECT_NAME} PRIVATE 
#     ${Python3_LIBRARIES}
#     ${FFTW3_LIBRARIES}
# )

# # Modern CMake target approach (if pkg-config succeeded)
# if(TARGET PkgConfig::FFTW3)
#     target_link_libraries(${PROJECT_NAME} PRIVATE PkgConfig::FFTW3)
# endif()

cmake_minimum_required(VERSION 3.24)
project(GPESolver VERSION 1.0.0 LANGUAGES CXX)

# —————————————————————————————————————————
# 1. C++ standard
# —————————————————————————————————————————
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# —————————————————————————————————————————
# 2. Platform-specific configuration
# —————————————————————————————————————————
if(APPLE)
  message(STATUS "Configuring for macOS (Homebrew)")
  
  # Python3 (Homebrew)
  set(Python3_ROOT_DIR "/Library/Frameworks/Python.framework/Versions/3.11")
  set(Python3_EXECUTABLE "${Python3_ROOT_DIR}/bin/python3.11")
  find_package(Python3 3.11 REQUIRED COMPONENTS Interpreter Development NumPy)
  
  # Eigen3 (Homebrew)
  find_path(EIGEN3_INCLUDE_DIR NAMES Eigen/Core
      PATHS /opt/homebrew/include /usr/local/include /usr/include
      PATH_SUFFIXES eigen3)
  if(NOT EIGEN3_INCLUDE_DIR)
      message(FATAL_ERROR "Eigen3 not found on macOS. Please brew install eigen")
  endif()
  
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
  message(STATUS "Configuring for PL-Grid (Linux)")
  
  # Python3
  find_package(Python3 3.11 REQUIRED COMPONENTS Interpreter Development NumPy)
  
  # Eigen3 (header-only)
  find_package(Eigen3 3.3 REQUIRED NO_MODULE)
  
  # FFTW3
  find_package(FFTW3 REQUIRED)
  
  # Threads (for fftw3_threads if needed)
  find_package(Threads)
endif()

# —————————————————————————————————————————
# 3. Executable and linking
# —————————————————————————————————————————
add_executable(${PROJECT_NAME} main.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE
  ${Python3_INCLUDE_DIRS}
  ${Python3_NumPy_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}            # for matplotlibcpp.h
  $<$<BOOL:${EIGEN3_INCLUDE_DIR}>:${EIGEN3_INCLUDE_DIR}>  # macOS Eigen or PL-Grid Eigen3::Eigen
  $<$<TARGET_EXISTS:Eigen3::Eigen>:${Eigen3_INCLUDE_DIR}> # PL-Grid Eigen
  $<$<BOOL:${FFTW3_INCLUDE_DIRS}>:${FFTW3_INCLUDE_DIRS}>  # macOS FFTW
  $<$<TARGET_EXISTS:FFTW3::fftw3>:${FFTW3_INCLUDE_DIRS}>  # PL-Grid FFTW
)

target_link_libraries(${PROJECT_NAME} PRIVATE
  Python3::Python
  Python3::NumPy
  $<$<TARGET_EXISTS:Eigen3::Eigen>:Eigen3::Eigen>
  $<$<TARGET_EXISTS:FFTW3::fftw3>:FFTW3::fftw3>
  $<$<BOOL:${FFTW3_LIBRARIES}>:${FFTW3_LIBRARIES}>  # macOS FFTW lib
#   Threads::Threads
)

# —————————————————————————————————————————
# 4. Final messages
# —————————————————————————————————————————
message(STATUS "Configuration complete. Build with: mkdir build && cd build && cmake .. && make -j")