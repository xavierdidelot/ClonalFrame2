cmake_minimum_required(VERSION 4.0)
project(cf2)

set(CMAKE_CXX_STANDARD 20)

add_library(cf2_lib
        Param.cpp
        Param.h
        Move.cpp
        Move.h)

add_executable(cf2 main.cpp)
target_link_libraries(cf2 PRIVATE cf2_lib)

