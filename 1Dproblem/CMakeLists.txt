cmake_minimum_required(VERSION 3.27.5)

set(PROJECT_NAME 1Dproblem)
project(${PROJECT_NAME} VERSION 1.0.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Python configuration
set(Python3_ROOT_DIR "/Library/Frameworks/Python.framework/Versions/3.11")
set(Python3_EXECUTABLE "/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11")

find_package(Python3 3.11 REQUIRED COMPONENTS Interpreter Development NumPy)

message(STATUS "Python executable: ${Python3_EXECUTABLE}")
message(STATUS "Python include directory: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python NumPy include directory: ${Python3_NumPy_INCLUDE_DIRS}")
message(STATUS "Python libraries: ${Python3_LIBRARIES}")

# Eigen configuration
find_path(EIGEN3_INCLUDE_DIR NAMES Eigen/Core
    PATHS /opt/homebrew/include /usr/include /usr/local/include
    PATH_SUFFIXES eigen3)

if(NOT EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR "Can't find Eigen")
endif()

message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}")

# # FFTW3 configuration (try pkg-config first, then fallback)
# find_package(PkgConfig REQUIRED)
# pkg_check_modules(FFTW3 IMPORTED_TARGET fftw3)

# if(FFTW3_FOUND)
#     message(STATUS "FFTW3 found via pkg-config")
#     message(STATUS "FFTW3 include dirs: ${FFTW3_INCLUDE_DIRS}")
#     message(STATUS "FFTW3 libraries: ${FFTW3_LIBRARIES}")
# else()
# message(WARNING "FFTW3 not found via pkg-config, trying manual search")

# Manual search for FFTW (common Homebrew locations)
find_path(FFTW3_INCLUDE_DIR fftw3.h
    PATHS /opt/homebrew/include /usr/local/include /usr/include
)

find_library(FFTW3_LIBRARY
    NAMES fftw3
    PATHS /opt/homebrew/lib /usr/local/lib /usr/lib
)

if(FFTW3_INCLUDE_DIR AND FFTW3_LIBRARY)
    message(STATUS "FFTW3 found manually")
    set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
    set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})
else()
    message(FATAL_ERROR "FFTW3 not found - please install via 'brew install fftw'")
endif()
# endif()

# Create executable
add_executable(${PROJECT_NAME} main.cpp ./CrankNicolson/CrankNicolson.cpp)

# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE
    ${Python3_INCLUDE_DIRS}        # Python 3.11
    ${Python3_NumPy_INCLUDE_DIRS}  # NumPy
    ${EIGEN3_INCLUDE_DIR}          # Eigen
    ${CMAKE_SOURCE_DIR}            # matplotlibcpp.h
    ${FFTW3_INCLUDE_DIRS}          # FFTW3
)

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE 
    ${Python3_LIBRARIES}
    ${FFTW3_LIBRARIES}
)

# Modern CMake target approach (if pkg-config succeeded)
if(TARGET PkgConfig::FFTW3)
    target_link_libraries(${PROJECT_NAME} PRIVATE PkgConfig::FFTW3)
endif()



# cmake_minimum_required(VERSION 3.27.5)

# set(PROJECT_NAME 1Dproblem)
# project(${PROJECT_NAME} VERSION 1.0.0)

# #set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED True)

# set(Python3_ROOT_DIR "/Library/Frameworks/Python.framework/Versions/3.11")
# set(Python3_EXECUTABLE "/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11")



# find_package(Python3 3.11 REQUIRED COMPONENTS Interpreter Development NumPy)

# message(STATUS "Python executable: ${Python3_EXECUTABLE}")
# message(STATUS "Python include directory: ${Python3_INCLUDE_DIRS}")
# message(STATUS "Python NumPy include directory: ${Python3_NumPy_INCLUDE_DIRS}")
# message(STATUS "Python libraries: ${Python3_LIBRARIES}")

# find_path(EIGEN3_INCLUDE_DIR NAMES Eigen/Core
#     PATHS /opt/homebrew/include /usr/include /usr/local/include
#     PATH_SUFFIXES eigen3)

# if(NOT EIGEN3_INCLUDE_DIR)
#     message(FATAL_ERROR "Can't find Eigen")
# endif()

# message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}")

# add_executable(${PROJECT_NAME} main.cpp ./CrankNicolson/CrankNicolson.cpp)

# target_include_directories(${PROJECT_NAME} PRIVATE
#     ${Python3_INCLUDE_DIRS}        # Python 3.11
#     ${Python3_NumPy_INCLUDE_DIRS}  # NumPy
#     ${EIGEN3_INCLUDE_DIR}          # Eigen
#     ${CMAKE_SOURCE_DIR}            # matplotlibcpp.h
# )
# target_link_libraries(${PROJECT_NAME} PRIVATE ${Python3_LIBRARIES})

