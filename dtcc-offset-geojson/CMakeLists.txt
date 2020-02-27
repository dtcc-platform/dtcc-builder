set(PROGRAM vc-offset-geojson)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM} pugixml proj shp las)
install(TARGETS ${PROGRAM} DESTINATION bin)
