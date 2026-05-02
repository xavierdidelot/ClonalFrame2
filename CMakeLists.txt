cmake_minimum_required(VERSION 3.20)
project(ClonalFrame2 VERSION 0.0.1 LANGUAGES CXX)
find_package(OpenMP REQUIRED)

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

add_subdirectory(src)
add_subdirectory(tests)

add_custom_target(doc
        COMMAND doxygen ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generation API documentation with Doxygen"
        VERBATIM
)