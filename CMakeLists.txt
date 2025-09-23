cmake_minimum_required(VERSION 4.0)
project(ClonalFrame2)
include(CTest)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_subdirectory(src)
enable_testing()
add_subdirectory(tests)
