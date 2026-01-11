add_executable(${PROJECT_NAME}
  main.cpp
  cli/cli.cpp
)

target_link_libraries(${PROJECT_NAME} PRIVATE
  CLI11::CLI11
)

if(GPES_ENABLE_PLOTTING)
  target_link_libraries(${PROJECT_NAME} PRIVATE
    gpes::plot
  )
else()
  target_link_libraries(${PROJECT_NAME} PRIVATE
    gpes::gpes
  )
endif()


if(GPES_WARNINGS)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endif()
