cmake_minimum_required(VERSION 3.20)
set(CMAKE_C_COMPILER /usr/local/bin/gcc-15)
set(CMAKE_CXX_COMPILER /usr/local/bin/g++-15)
project(ClonalFrame2 VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(VERSION_STRING "${PROJECT_VERSION}")

configure_file(src/version.h.in ${CMAKE_BINARY_DIR}/generated/version.h)

include(CTest)
enable_testing()

add_library(cf2_lib
        src/Alignment.cpp
        src/Alignment.h
        src/Tree.cpp
        src/Tree.h
        src/Node.cpp
        src/Node.h
        src/Param.cpp
        src/Param.h
        src/Profiler.cpp
        src/Profiler.h)
add_executable(ClonalFrame2 src/main.cpp)
target_link_libraries(ClonalFrame2 PRIVATE cf2_lib)
target_include_directories(ClonalFrame2 PRIVATE ${CMAKE_BINARY_DIR}/generated)

add_executable(Alignment_test tests/Alignment_test.cpp)
target_link_libraries(Alignment_test PRIVATE cf2_lib)
add_test(Alignment_test Alignment_test)
add_executable(Tree_test tests/Tree_test.cpp)
target_link_libraries(Tree_test PRIVATE cf2_lib)
add_test(Tree_test Tree_test)
add_executable(Param_test tests/Param_test.cpp)
target_link_libraries(Param_test PRIVATE cf2_lib)
add_test(Param_test Param_test)
add_executable(test_nj tests/test_nj.cpp)
target_link_libraries(test_nj PRIVATE cf2_lib)
add_executable(test_parsimony tests/test_parsimony.cpp)
target_link_libraries(test_parsimony PRIVATE cf2_lib)
add_executable(test_patterns tests/test_patterns.cpp)
target_link_libraries(test_patterns PRIVATE cf2_lib)
add_executable(test_likelihoods tests/test_likelihoods.cpp)
target_link_libraries(test_likelihoods PRIVATE cf2_lib)
add_executable(test_spr tests/test_spr.cpp)
target_link_libraries(test_spr PRIVATE cf2_lib)
add_executable(test_jointml tests/test_jointml.cpp)
target_link_libraries(test_jointml PRIVATE cf2_lib)
add_executable(profile_parsimony tests/profile_parsimony.cpp)
target_link_libraries(profile_parsimony PRIVATE cf2_lib)
add_executable(comparison tests/comparison.cpp)
target_link_libraries(comparison PRIVATE cf2_lib)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(cf2_lib PUBLIC OpenMP::OpenMP_CXX)
else()
    message(WARNING "OpenMP not found - building without parallelisation")
endif()

add_custom_target(doc
        COMMAND doxygen ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generation API documentation with Doxygen"
        VERBATIM
)