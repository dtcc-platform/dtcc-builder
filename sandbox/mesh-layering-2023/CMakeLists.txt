cmake_minimum_required(VERSION 3.12)
project(mesh_layering)

find_package(pybind11 REQUIRED)

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cpp/)
set(BUILDER_DIR $ENV{HOME}/dtcc-builder)

# Set include directories
include_directories(${BUILDER_DIR}/src/cpp/include/)
include_directories(${BUILDER_DIR}/src/cpp/include/model/)
include_directories(${BUILDER_DIR}/src/cpp/external)
include_directories(/usr/local/include)

pybind11_add_module(_mesh_layering ${SRC_DIR}/mesh_layering.cpp)