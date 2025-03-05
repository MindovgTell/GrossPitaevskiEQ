# cmake_minimum_required(VERSION 3.27.5)

# set(PROJECT_NAME 1Dproblem)

# project(${PROJECT_NAME} VERSION 1.0.0)

# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED True)

# # Find Eigen
# include_directories("/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")

# add_executable(${PROJECT_NAME} main.cpp ./CrankNicolson/CrankNicolson.cpp)
# target_include_directories(${PROJECT_NAME} PRIVATE "/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")




#________________________________________________________________________________________
# # cmake_minimum_required(VERSION 3.27.5)

# set(PROJECT_NAME 1Dproblem)
# project(${PROJECT_NAME} VERSION 1.0.0)

# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED True)

# # Find Eigen
# include_directories("/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")

# # Найти Python и NumPy
# find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

# # Подключаем директорию NumPy
# include_directories("/opt/homebrew/lib/python3.12/site-packages/numpy/_core/include")  # Путь к NumPy

# # Подключаем matplotlibcpp
# include_directories("./matplotlib-cpp")  # Путь к matplotlibcpp.h

# # Ваши исходные файлы
# add_executable(${PROJECT_NAME} main.cpp ./CrankNicolson/CrankNicolson.cpp)

# # Линковка с Python
# target_include_directories(${PROJECT_NAME} PRIVATE ${Python3_INCLUDE_DIRS})
# target_link_libraries(${PROJECT_NAME} PRIVATE ${Python3_LIBRARIES})

# # Подключаем директории для Eigen и matplotlibcpp
# target_include_directories(${PROJECT_NAME} PRIVATE "/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")
#________________________________________________________________________________________


cmake_minimum_required(VERSION 3.27.5)

# Название проекта
set(PROJECT_NAME 1Dproblem)
project(${PROJECT_NAME} VERSION 1.0.0)

# Стандарт C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Указываем путь к Python 3.11
set(Python3_EXECUTABLE "/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11")  # Укажите путь к Python 3.11

# Находим компоненты Python: интерпретатор, разработки и библиотеки
find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)

# Выводим путь к Python для проверки
message(STATUS "Python executable: ${Python3_EXECUTABLE}")
message(STATUS "Python include directory: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python libraries: ${Python3_LIBRARIES}")

# Устанавливаем путь к Eigen
include_directories("/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")

# Подключаем пути к библиотекам Python и NumPy
include_directories(${Python3_INCLUDE_DIRS})
include_directories("/Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/numpy/_core/include")

# Подключаем matplotlibcpp (укажите правильный путь к библиотеке)
include_directories("./")  # Укажите путь к matplotlibcpp.h

# Добавляем исходные файлы проекта
add_executable(${PROJECT_NAME} main.cpp ./CrankNicolson/CrankNicolson.cpp)

# Линковка с библиотеками Python
target_link_libraries(${PROJECT_NAME} PRIVATE ${Python3_LIBRARIES})

# Линковка с дополнительными библиотеками, если нужно (например, Eigen)
target_include_directories(${PROJECT_NAME} PRIVATE "/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")

