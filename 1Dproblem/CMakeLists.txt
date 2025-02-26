cmake_minimum_required(VERSION 3.27.5)

set(PROJECT_NAME 1Dproblem)

project(${PROJECT_NAME} VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find Eigen
include_directories("/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")

add_executable(${PROJECT_NAME} main.cpp ./CrankNicolson/CrankNicolson.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE "/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")

