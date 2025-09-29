add_library(cf2_lib
        Param.cpp
        Param.h
        Move.cpp
        Move.h
        Alignment.cpp
        Alignment.h
        Tree.cpp
        Tree.h
        Node.cpp
        Node.h)

add_executable(ClonalFrame2 main.cpp)
find_package(GSL REQUIRED)
target_link_libraries(ClonalFrame2 cf2_lib GSL::gsl GSL::gslcblas)

