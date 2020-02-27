set(PROGRAM dtcc-read-cdf)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM} pugixml proj shp las netcdf_c++4 netcdf)
install(TARGETS ${PROGRAM} DESTINATION bin)
