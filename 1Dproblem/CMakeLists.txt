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
# cmake_minimum_required(VERSION 3.27.5)

# set(PROJECT_NAME 1Dproblem)
# project(${PROJECT_NAME} VERSION 1.0.0)

# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED True)

# # Find Eigen
# include_directories("/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")

# # Найти Python и NumPy
# find_package(Python3.11 COMPONENTS Interpreter Development NumPy REQUIRED)

# # Подключаем директорию NumPy
# include_directories("/Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/numpy/core/include")  # Путь к NumPy
# include_directories("/opt/homebrew/Cellar/python-matplotlib/3.10.1")

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

# Указываем стандарт C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Жестко указываем путь к Python 3.11
set(Python3_ROOT_DIR "/Library/Frameworks/Python.framework/Versions/3.11")
set(Python3_EXECUTABLE "/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11")

# Находим Python 3.11 и NumPy
find_package(Python3 3.11 REQUIRED COMPONENTS Interpreter Development NumPy)

# Вывод информации о найденных компонентах Python
message(STATUS "Python executable: ${Python3_EXECUTABLE}")
message(STATUS "Python include directory: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python NumPy include directory: ${Python3_NumPy_INCLUDE_DIRS}")
message(STATUS "Python libraries: ${Python3_LIBRARIES}")

# Подключаем Eigen (Homebrew путь, проверь свою версию)
find_path(EIGEN3_INCLUDE_DIR NAMES Eigen/Core
    PATHS /opt/homebrew/include /usr/include /usr/local/include
    PATH_SUFFIXES eigen3)

if(NOT EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR "Eigen3 не найден. Установите его через Homebrew: brew install eigen")
endif()

message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}")

# Создаем исполняемый файл
add_executable(${PROJECT_NAME} main.cpp ./CrankNicolson/CrankNicolson.cpp)

# Добавляем заголовочные файлы
target_include_directories(${PROJECT_NAME} PRIVATE
    ${Python3_INCLUDE_DIRS}        # Python 3.11 заголовки
    ${Python3_NumPy_INCLUDE_DIRS}  # NumPy заголовки
    ${EIGEN3_INCLUDE_DIR}          # Eigen
    ${CMAKE_SOURCE_DIR}            # matplotlibcpp.h (если он в корне проекта)
)

# Линковка с Python 3.11
target_link_libraries(${PROJECT_NAME} PRIVATE ${Python3_LIBRARIES})


