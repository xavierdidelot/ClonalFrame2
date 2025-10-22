##find_package(GSL REQUIRED)

add_executable(Alignment_test Alignment_test.cpp)
target_link_libraries(Alignment_test cf2_lib)# GSL::gsl GSL::gslcblas)
add_test(Alignment_test Alignment_test)

add_executable(Tree_test Tree_test.cpp)
target_link_libraries(Tree_test cf2_lib)# GSL::gsl GSL::gslcblas)
add_test(Tree_test Tree_test)

add_executable(Param_test Param_test.cpp)
target_link_libraries(Param_test cf2_lib)# GSL::gsl GSL::gslcblas)
add_test(Param_test Param_test)
