# Project configuration
cmake_minimum_required(VERSION 3.5)
project(Core)

# Note: To building in Debug mode, use cmake .. -DCMAKE_BUILD_TYPE=Debug

# Set default build type to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  message(STATUS "Build type not specified, using Release by default")
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

# Set build flags
set(CMAKE_CXX_FLAGS_ "-Wall -Werror")
set(CMAKE_CXX_FLAGS_DEBUG "-pg -DLLVM_ENABLE_ASSERTIONS=ON")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

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
find_package(nlohmann_json REQUIRED)
#find_package(assimp REQUIRED)
find_package(GDAL REQUIRED)
find_package(PNG REQUIRED)

# Find FEniCS
find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

# Hack to remove -DNDEBUG from FEniCS flags (shouldn't be there!)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  get_target_property(_old_flags dolfin INTERFACE_COMPILE_DEFINITIONS)
  set(_new_flags "")
  foreach (_old_flag ${_old_flags})
    if (NOT "${_old_flag}" STREQUAL "NDEBUG")
      list(APPEND _new_flags ${_old_flag})
    endif()
  endforeach()
  set_target_properties(dolfin PROPERTIES INTERFACE_COMPILE_DEFINITIONS "${_new_flags}")
endif()

# Find VTK
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# Set include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include/citymodel)
include_directories(${CMAKE_SOURCE_DIR}/external)
include_directories(${CMAKE_SOURCE_DIR}/external)
include_directories(/usr/include/nlohmann/)

link_directories(/usr/local/lib) #Equivalent of ldconfig /usr/local/lib in Dockerfile

# Set installation directory
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)

# Add subdirectories
add_subdirectory(dtcc-info)
add_subdirectory(dtcc-plot)
add_subdirectory(dtcc-generate-citymodel)
add_subdirectory(dtcc-generate-mesh)
add_subdirectory(dtcc-generate-roadnetwork)
add_subdirectory(dtcc-randomize-citymodel)
add_subdirectory(dtcc-smooth-ground-mesh)
add_subdirectory(dtcc-color-pointcloud)
#add_subdirectory(dtcc-convert-shp2geojson)
#add_subdirectory(dtcc-convert-geojson2shp)
add_subdirectory(dtcc-assimp-hello-world)
add_subdirectory(dtcc-assimp-normals)
add_subdirectory(dtcc-generate-cityjson-from-mesh)
#add_subdirectory(dtcc-offset-geojson)
add_subdirectory(dtcc-read-netcdf)
add_subdirectory(dtcc-sandbox)
#add_subdirectory(dtcc-iotools)
add_subdirectory(unittests)

# Print all variables (uncomment for CMake debugging)
#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

# Pre-commit hook
MESSAGE("Running pre-commit hook")
execute_process(COMMAND bash -c "${CMAKE_SOURCE_DIR}/scripts/git/git-pre-commit-format install"
                OUTPUT_VARIABLE outVar)
MESSAGE("${outVar}")
