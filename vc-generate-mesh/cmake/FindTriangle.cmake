# FindTriangle.cmake (by Anders Logg 2018)
#
# Try to find Triangle and set the following variables:
#
#   Triangle_FOUND
#   Triangle_INCLUDE_DIRS
#   Triangle_LIBRARIES

find_path(Triangle_INCLUDE_DIR
  NAMES triangle.h
  PATHS /usr/local /usr
  PATH_SUFFIXES include
)

find_library(Triangle_LIBRARY
  NAMES triangle
  PATHS /usr/local /usr
  PATH_SUFFIXES lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Triangle
  FOUND_VAR Triangle_FOUND
  REQUIRED_VARS
    Triangle_INCLUDE_DIR
    Triangle_LIBRARY
)

if(Triangle_FOUND)
  set(Triangle_INCLUDE_DIRS ${Triangle_INCLUDE_DIR})
  set(Triangle_LIBRARIES ${Triangle_LIBRARY})
endif()

mark_as_advanced(
  Triangle_INCLUDE_DIR
  Triangle_LIBRARY
)
