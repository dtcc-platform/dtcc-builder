# Find the PROJ4 library and set the following standard variables:
#
# PROJ4_INCLUDE_DIRS - header location
# PROJ4_LIBRARIES    - library to link against
# PROJ4_FOUND        - true if library was found

include(FindPackageHandleStandardArgs)

find_path(PROJ4_INCLUDE_DIRS
          NAMES proj_api.h
          PATHS /usr/local/include)

find_library(PROJ4_LIBRARIES
             NAMES proj
             PATHS /usr/local/lib)

find_package_handle_standard_args(PROJ4
                                  DEFAULT_MSG
                                  PROJ4_INCLUDE_DIRS
                                  PROJ4_LIBRARIES)
