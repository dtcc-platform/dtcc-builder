include(FindPackageHandleStandardArgs)

find_path(DTCCMODEL_INCLUDE_DIRS
          NAMES dtcc.pb.h
          PATHS /usr/local/include/dtcc /usr/local/include /usr/include/dtcc /usr/include)

find_library(DTCCMODEL_LIBRARIES
             NAMES dtccproto
             PATHS /usr/local/lib /usr/lib)

find_package_handle_standard_args(DTCCModel
                                  DEFAULT_MSG
                                  DTCCMODEL_INCLUDE_DIRS
                                  DTCCMODEL_LIBRARIES)