find_package(GSL REQUIRED)

add_executable(Param_test Param_test.cpp)
target_link_libraries(Param_test cf2_lib GSL::gsl GSL::gslcblas)
add_test(Param_test Param_test)

add_executable(Alignment_test Alignment_test.cpp)
target_link_libraries(Alignment_test cf2_lib GSL::gsl GSL::gslcblas)
add_test(Alignment_test Alignment_test)
