add_library(cf2_lib
        Alignment.cpp
        Alignment.h
        Tree.cpp
        Tree.h
        Node.cpp
        Node.h
        Param.cpp
        Param.h)

add_executable(ClonalFrame2 main.cpp)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(ClonalFrame2 OpenMP::OpenMP_CXX)
endif()
target_link_libraries(ClonalFrame2 cf2_lib)
target_include_directories(ClonalFrame2 PRIVATE ${CMAKE_BINARY_DIR}/generated)
