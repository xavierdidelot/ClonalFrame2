cmake_minimum_required(VERSION 4.0)
project(cf2)

set(CMAKE_CXX_STANDARD 20)

add_executable(cf2 main.cpp
        Param.cpp
        Param.h
        Move.cpp
        Move.h)
