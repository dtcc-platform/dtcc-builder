add_executable(vc-convert-shp2geojson vc-convert-shp2geojson.cpp)
target_link_libraries(vc-convert-shp2geojson pugixml proj shp las)
install(TARGETS vc-convert-shp2geojson DESTINATION bin)
