add_executable(Param_test Param_test.cpp)
target_link_libraries(Param_test PRIVATE cf2_lib)
add_test(Param_test Param_test)
