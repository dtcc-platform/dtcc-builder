# Project configuration
cmake_minimum_required(VERSION 3.12)
project(dtcc-builder)

# Note: To building in Debug mode, use cmake .. -DCMAKE_BUILD_TYPE=Debug
# Set default build type to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  message(STATUS "Build type not specified, using Release by default")
endif()

# MacOS doesn't ship with OpenMP, so we need to use Homebrew's Clang
# if we want OpenMP support. To get OpenMP support on Mac run:
# brew install llvm
# brew install libomp
if(APPLE)
  find_program(HOMEBREW_CLANG /opt/homebrew/opt/llvm/bin/clang)
  find_program(HOMEBREW_CLANGXX /opt/homebrew/opt/llvm/bin/clang++)

  # Set Homebrew compiler paths, if found, and set OpenMP flags
  if(HOMEBREW_CLANG AND HOMEBREW_CLANGXX)
    # set(CMAKE_PREFIX_PATH $(brew --prefix libomp))
    # Get Homebrew prefix path for libomp
    execute_process(COMMAND brew --prefix libomp OUTPUT_VARIABLE HOMEBREW_LIBOMP_PREFIX)
    string(STRIP ${HOMEBREW_LIBOMP_PREFIX} HOMEBREW_LIBOMP_PREFIX)
    list(APPEND CMAKE_PREFIX_PATH ${HOMEBREW_LIBOMP_PREFIX})
    message(STATUS "Using Homebrew Clang (${HOMEBREW_CLANG}) and Clang++ (${HOMEBREW_CLANGXX}) compilers")
    message(STATUS "CMake prefix ${CMAKE_PREFIX_PATH}")
    set(CMAKE_C_COMPILER ${HOMEBREW_CLANG})
    set(CMAKE_CXX_COMPILER ${HOMEBREW_CLANGXX})
  else()
    message(WARNING "Homebrew Clang and Clang++ compilers not found. Using default compilers. 
                      OpenMP will be disabled.")
  endif()
endif()

message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

# Set build flags
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-reorder")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg -DLLVM_ENABLE_ASSERTIONS=ON")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
endif()

# Add path to local CMake scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Add clang tools
include(cmake/clang-cxx-dev-tools.cmake)

# # Find pybind11
execute_process(COMMAND python -c "import site; print([s for s in site.getsitepackages() if 'site-packages' in s][0])" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
set(pybind11_DIR ${PYTHON_SITE_PACKAGES}/pybind11/share/cmake/pybind11)

# set(PYBIND11_FINDPYTHON ON)
message("Pybind11 DIR " ${pybind11_DIR})
find_package(pybind11 REQUIRED)

# Set include directories
include_directories(${CMAKE_SOURCE_DIR}/src/cpp/include)
include_directories(${CMAKE_SOURCE_DIR}/src/cpp/include/citymodel)
include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/src/cpp/external)
include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/src/cpp/external/triangle)

# Set installation directory
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/src)
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)

# Add subdirectories
# add_subdirectory(src/cpp/external/triangle)
add_subdirectory(src)

if(NOT WIN32)
  # tests won't build on Windows
  # TODO: Work out why
  add_subdirectory(tests/cpp)
endif()

# Pre-commit hook
MESSAGE("Running pre-commit hook")
execute_process(COMMAND bash -c "${CMAKE_SOURCE_DIR}/utils/git/git-pre-commit-format install"
  OUTPUT_VARIABLE outVar)
MESSAGE("${outVar}")
