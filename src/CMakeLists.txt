cmake_minimum_required (VERSION 2.8)

project (HELWFEM)

set(CMAKE_BUILD_TYPE debug)
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++0x ${CMAKE_CXX_FLAGS} -Wextra -g")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++0x ${CMAKE_CXX_FLAGS} -Wextra -g -O3")

set(CMAKE_BINARY_DIR build)
set(EXECUTABLE_OUTPUT_PATH bin)

add_executable (helwfem main.cpp assembly.cpp element.cpp boundary.cpp)

target_link_libraries (helwfem armadillo)

install (TARGETS helwfem DESTINATION bin)

