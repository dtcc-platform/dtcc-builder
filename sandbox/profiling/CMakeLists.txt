# Project configuration
cmake_minimum_required(VERSION 3.5)
project(dtcc-builder VERSION 0.1)

# Set the output directory for the generated Makefile
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Set build flags
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-std=c++14 -Wall -Wextra -Wpedantic -O3")
# set(CMAKE_CXX_FLAGS_DEBUG "-Werror -pg -DLLVM_ENABLE_ASSERTIONS=ON")
# set(CMAKE_CXX_FLAGS_RELEASE " -DNDEBUG -Wunused-function -Wunused-variable")


# Add clang tools
# include(cmake/clang-cxx-dev-tools.cmake)
# include(cmake/FindLibGEOS.cmake)

# Using Clang-tidy for static analysis
set(CMAKE_CXX_CLANG_TIDY "clang-tidy;-checks=clang-analyzer-* -extra-arg=-std=c++14  -extra-arg=-lstdc++fs --fix-errors ")

# Find packages


# Set include directories
include_directories(/home/auth-georspai/dtcc-builder/src/cpp/include/)
include_directories(/home/auth-georspai/dtcc-builder/src/cpp/include/model/)
include_directories(/home/auth-georspai/dtcc-builder/src/cpp/external)
include_directories(/usr/local/include)

# find_package(Triangle REQUIRED)
# find_package(LibGEOS 3.16 REQUIRED)
# Set installation directory
# set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
# set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)

set(SOURCES SmootherProfiling.cpp)

add_executable(smoother-profiling ${SOURCES})

