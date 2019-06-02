add_executable(vc-offset-geojson vc-offset-geojson.cpp)
target_link_libraries(vc-offset-geojson pugixml proj shp las)
install(TARGETS vc-offset-geojson DESTINATION bin)
