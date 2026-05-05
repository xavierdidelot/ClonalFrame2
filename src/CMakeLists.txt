add_library(cf2_lib
        Alignment.cpp
        Alignment.h
        Tree.cpp
        Tree.h
        Node.cpp
        Node.h
        Param.cpp
        Param.h
        Profiler.cpp
        Profiler.h)
target_compile_options(cf2_lib PRIVATE -fopenmp)

add_executable(ClonalFrame2 main.cpp)
target_link_libraries(ClonalFrame2 PRIVATE OpenMP::OpenMP_CXX cf2_lib)
target_include_directories(ClonalFrame2 PRIVATE ${CMAKE_BINARY_DIR}/generated)
