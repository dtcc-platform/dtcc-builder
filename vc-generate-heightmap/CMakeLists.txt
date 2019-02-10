cmake_minimum_required(VERSION 3.5)
project(VCCore)
set(CMAKE_BUILD_TYPE Debug)

# Add path to local CMake scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Add path to common includes
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

# Find ImageMagick
find_package(ImageMagick COMPONENTS Magick++ REQUIRED)
include_directories(${ImageMagick_INCLUDE_DIRS})

# Build executable
add_executable(vc-generate-heightmap vc-generate-heightmap.cpp)
target_link_libraries(vc-generate-heightmap ${ImageMagick_LIBRARIES})
