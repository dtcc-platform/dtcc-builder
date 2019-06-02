add_executable(vc-convert-geojson2shp vc-convert-geojson2shp.cpp)
target_link_libraries(vc-convert-geojson2shp pugixml proj shp las)
install(TARGETS vc-convert-geojson2shp DESTINATION bin)
