# Find the LibLAS library and set the following standard variables:
#
# LIBLAS_INCLUDE_DIRS - header location
# LIBLAS_LIBRARIES    - library to link against
# LIBLAS_FOUND        - true if library was found

include(FindPackageHandleStandardArgs)

find_path(LIBLAS_INCLUDE_DIRS
          NAMES liblas/liblas.hpp
          PATHS /usr/local/include)

find_library(LIBLAS_LIBRARIES
             NAMES las
             PATHS /usr/local/lib)

find_package_handle_standard_args(LIBLAS
                                  DEFAULT_MSG
                                  LIBLAS_INCLUDE_DIRS
                                  LIBLAS_LIBRARIES)

