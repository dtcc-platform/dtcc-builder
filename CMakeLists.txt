# Project configuration
cmake_minimum_required(VERSION 3.12)
project(dtcc-builder)

# Note: To building in Debug mode, use cmake .. -DCMAKE_BUILD_TYPE=Debug

# Set default build type to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  message(STATUS "Build type not specified, using Release by default")
endif()

message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

# Set build flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_ "-Wall -Werror")
set(CMAKE_CXX_FLAGS_DEBUG "-pg -DLLVM_ENABLE_ASSERTIONS=ON")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")


# Add path to local CMake scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Add clang tools
include(cmake/clang-cxx-dev-tools.cmake)

## Find pybind11

execute_process (COMMAND python -c "import site; print(site.getsitepackages()[0])" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
set(pybind11_DIR ${PYTHON_SITE_PACKAGES}/pybind11/share/cmake/pybind11)
set(PYBIND11_FINDPYTHON ON)

message("Pybind11 DIR " ${pybind11_DIR})
find_package(pybind11 REQUIRED)

# Set include directories
include_directories(${CMAKE_SOURCE_DIR}/src/cpp/include)
include_directories(${CMAKE_SOURCE_DIR}/src/cpp/include/citymodel)
include_directories(${CMAKE_SOURCE_DIR}/src/cpp/external)
include_directories(${CMAKE_SOURCE_DIR}/src/cpp/external/triangle)

# Set installation directory
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/src)
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)


# Add subdirectories
# add_subdirectory(src/cpp/external/triangle)
add_subdirectory(src)

if (NOT WIN32)
  # tests won't build on Windows
  # TODO: Work out why
  add_subdirectory(tests/cpp)
endif()

# Pre-commit hook
MESSAGE("Running pre-commit hook")
execute_process(COMMAND bash -c "${CMAKE_SOURCE_DIR}/utils/git/git-pre-commit-format install"
  OUTPUT_VARIABLE outVar)
MESSAGE("${outVar}")
