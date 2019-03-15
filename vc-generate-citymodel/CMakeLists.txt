cmake_minimum_required(VERSION 3.5)
project(VCCore)
set(CMAKE_BUILD_TYPE Debug)

# Add path to local CMake scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Add path to common includes
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

# Find packages
find_package(PugiXML REQUIRED)
find_package(PROJ4 REQUIRED)
find_package(ShapeLib REQUIRED)

# Build executable
add_executable(vc-generate-citymodel vc-generate-citymodel.cpp)
target_link_libraries(vc-generate-citymodel pugixml proj shp)
