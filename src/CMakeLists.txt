cmake_minimum_required (VERSION 2.8)

project (HELWFEM)

SET(CMAKE_BUILD_TYPE debug)
set(CMAKE_C_FLAGS_DEBUG "-Wextra -g")
set(CMAKE_C_FLAGS_RELEASE "-Wextra -g -O3")

set(CMAKE_BINARY_DIR build)
set(EXECUTABLE_OUTPUT_PATH bin)

add_executable (helwfem main.cpp assembly.cpp element.cpp)

target_link_libraries (helwfem armadillo)

install (TARGETS helwfem DESTINATION bin)

