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
#find_package(GSL REQUIRED)
target_link_libraries(ClonalFrame2 cf2_lib)# GSL::gsl GSL::gslcblas)
target_include_directories(ClonalFrame2 PRIVATE ${CMAKE_BINARY_DIR}/generated)
