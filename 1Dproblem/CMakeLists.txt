cmake_minimum_required(VERSION 3.27.5)

set(PROJECT_NAME 1Dproblem)
project(${PROJECT_NAME} VERSION 1.0.0)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(Python3_ROOT_DIR "/Library/Frameworks/Python.framework/Versions/3.11")
set(Python3_EXECUTABLE "/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11")

find_package(Python3 3.11 REQUIRED COMPONENTS Interpreter Development NumPy)

message(STATUS "Python executable: ${Python3_EXECUTABLE}")
message(STATUS "Python include directory: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python NumPy include directory: ${Python3_NumPy_INCLUDE_DIRS}")
message(STATUS "Python libraries: ${Python3_LIBRARIES}")

find_path(EIGEN3_INCLUDE_DIR NAMES Eigen/Core
    PATHS /opt/homebrew/include /usr/include /usr/local/include
    PATH_SUFFIXES eigen3)

if(NOT EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR "Can't find Eigen")
endif()

message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}")

add_executable(${PROJECT_NAME} main.cpp ./CrankNicolson/CrankNicolson.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE
    ${Python3_INCLUDE_DIRS}        # Python 3.11
    ${Python3_NumPy_INCLUDE_DIRS}  # NumPy
    ${EIGEN3_INCLUDE_DIR}          # Eigen
    ${CMAKE_SOURCE_DIR}            # matplotlibcpp.h
)
target_link_libraries(${PROJECT_NAME} PRIVATE ${Python3_LIBRARIES})


