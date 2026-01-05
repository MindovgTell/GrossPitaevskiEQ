add_executable(${PROJECT_NAME}
  main.cpp
)

target_link_libraries(${PROJECT_NAME} PRIVATE
  gpes::gpes
  CLI11::CLI11
)

# # matplotlib-cpp (header-only), if you use it in the CLI
# if(GPES_ENABLE_PYTHON_PLOTTING)
#   target_link_libraries(${PROJECT_NAME} PRIVATE
#     Python3::Python
#     Python3::NumPy
#   )

#   target_include_directories(${PROJECT_NAME} PRIVATE
#     ${Python3_INCLUDE_DIRS}
#     ${Python3_NumPy_INCLUDE_DIRS}
#     ${CMAKE_SOURCE_DIR}/third_party/matplotlib-cpp
#   )
# endif()

if(GPES_WARNINGS)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endif()