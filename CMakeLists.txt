# Project configuration
cmake_minimum_required(VERSION 3.5)
project(Core)
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
find_package(netCDF REQUIRED)
find_package(VTK REQUIRED)
find_package(assimp REQUIRED)
# Find FEniCS
find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

# Set include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/external)

# Set installation directory
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

# Add subdirectories
add_subdirectory(dtcc-info)
add_subdirectory(dtcc-plot)
add_subdirectory(dtcc-generate-heightmap)
add_subdirectory(dtcc-generate-citymodel)
add_subdirectory(dtcc-simplify-citymodel)
add_subdirectory(dtcc-generate-simulation-mesh)
add_subdirectory(dtcc-generate-visualization-mesh)
add_subdirectory(dtcc-smooth-ground-mesh)
#add_subdirectory(dtcc-convert-shp2geojson)
#add_subdirectory(dtcc-convert-geojson2shp)
add_subdirectory(dtcc-assimp-hello-world)
add_subdirectory(dtcc-assimp-normals)
#add_subdirectory(dtcc-offset-geojson)
add_subdirectory(dtcc-read-netcdf)
add_subdirectory(dtcc-sandbox)
#add_subdirectory(dtcc-iotools)
add_subdirectory(unittests) 
MESSAGE("Running pre-commit hook") 
execute_process(COMMAND bash -c "${CMAKE_SOURCE_DIR}/scripts/git-pre-commit-format install"
                OUTPUT_VARIABLE outVar) 
MESSAGE("${outVar}")
