# Project configuration
cmake_minimum_required(VERSION 3.5)
project(VCCore)
set(CMAKE_BUILD_TYPE Debug)

# Set build flags
set (CMAKE_CXX_FLAGS "-Wall -Werror")
#-Wno-sign-compare")

# Add path to local CMake scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Add clang tools
include(cmake/clang-cxx-dev-tools.cmake)

# Find packages
find_package(ShapeLib REQUIRED)
find_package(LibLAS)
find_package(PugiXML REQUIRED)
find_package(PROJ4 REQUIRED)
find_package(Triangle REQUIRED)
find_package(CDF REQUIRED)

# Find FEniCS
find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

# Set include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Set installation directory
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

# Add subdirectories
add_subdirectory(vc-info)
add_subdirectory(vc-plot)
add_subdirectory(vc-generate-heightmap)
add_subdirectory(vc-generate-citymodel)
add_subdirectory(vc-simplify-citymodel)
add_subdirectory(vc-generate-simulation-mesh)
add_subdirectory(vc-generate-visualization-mesh)
add_subdirectory(vc-convert-shp2geojson)
add_subdirectory(vc-convert-geojson2shp)
#add_subdirectory(vc-offset-geojson)
add_subdirectory(vc-read-netcdf)
