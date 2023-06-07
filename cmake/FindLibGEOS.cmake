# Find the LibLAS library and set the following standard variables:
#
# LIBGEOS_INCLUDE_DIRS - header location
# LIBGEOS_LIBRARIES    - library to link against
# LIBGEOS_FOUND        - true if library was found

include(FindPackageHandleStandardArgs)

find_path(LIBGEOS_INCLUDE_DIRS
    NAMES geos_c.h
    PATHS /usr/local/include /opt/homebrew/include)

find_library(LIBGEOS_LIBRARIES
    NAMES geos_c
    PATHS /usr/local/lib /opt/homebrew/lib)

find_package_handle_standard_args(LibGEOS
    DEFAULT_MSG
    LIBGEOS_INCLUDE_DIRS
    LIBGEOS_LIBRARIES)
