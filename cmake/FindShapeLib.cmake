# Find the ShapeLib library and set the following standard variables:
#
# SHAPELIB_INCLUDE_DIRS - header location
# SHAPELIB_LIBRARIES    - library to link against
# SHAPELIB_FOUND        - true if library was found

include(FindPackageHandleStandardArgs)

find_path(SHAPELIB_INCLUDE_DIRS
          NAMES shapefil.h
          PATHS /usr/local/include)

find_library(SHAPELIB_LIBRARIES
             NAMES shp
             PATHS /usr/local/lib)

find_package_handle_standard_args(SHAPELIB
                                  DEFAULT_MSG
                                  SHAPELIB_INCLUDE_DIRS
                                  SHAPELIB_LIBRARIES)
