add_executable(gpes_numerics_test numerics.cpp)
target_link_libraries(gpes_numerics_test PRIVATE gpes::gpes)
target_compile_options(gpes_numerics_test PRIVATE
  $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wall;-Wextra;-Wpedantic>
)
add_test(NAME gpes_numerics COMMAND gpes_numerics_test)
