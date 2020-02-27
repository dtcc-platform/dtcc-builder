cmake_minimum_required(VERSION 3.5)
project(Core)
set(CMAKE_BUILD_TYPE Debug)

# Add path to local CMake scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Add path to common includes
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

# Build executable
add_executable(dtcc-randomize-citymodel dtcc-randomize-citymodel.cpp)
target_link_libraries(dtcc-randomize-citymodel)
