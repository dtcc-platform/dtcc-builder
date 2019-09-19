# Find the PROJ4 library and set the following standard variables:
#
# PROJ4_INCLUDE_DIRS - header location
# PROJ4_LIBRARIES    - library to link against
# PROJ4_FOUND        - true if library was found

include(FindPackageHandleStandardArgs)

find_path(CDF_INCLUDE_DIRS
          NAMES netcdf.h
          PATHS /usr/local/include)

find_library(CDF_LIBRARIES
             NAMES netcdf_c++4 netcdf
             PATHS /usr/lib/x86_64-linux-gnu)
find_package_handle_standard_args(CDF
                                  DEFAULT_MSG
                                  CDF_INCLUDE_DIRS
                                  CDF_LIBRARIES)
