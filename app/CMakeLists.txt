add_executable(${PROJECT_NAME}
  main.cpp
  cli/cli.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../include/gpes
)

# Add executable-only compile flags here.
# Example:
#   -O3
#   -march=native
set(GPES_EXECUTABLE_COMPILE_OPTIONS
    -O3
    -march=native
    -ffast-math
    -s
)

target_compile_options(${PROJECT_NAME} PRIVATE
  ${GPES_EXECUTABLE_COMPILE_OPTIONS}
)

if(GPES_WARNINGS)
  target_compile_options(${PROJECT_NAME} PRIVATE
    $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wall;-Wextra;-Wpedantic>
  )
endif()

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
