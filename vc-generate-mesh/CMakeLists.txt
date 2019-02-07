cmake_minimum_required(VERSION 3.5)
project(vc-backend)
set(CMAKE_BUILD_TYPE Debug)

# Add path to local CMake scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Find DOLFIN
find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

# Find Triangle
find_package(Triangle REQUIRED)

# Build executable
add_executable(vc-generate-mesh vc-generate-mesh.cpp)
target_link_libraries(vc-generate-mesh dolfin triangle)
