add_executable(vc-generate-citymodel vc-generate-citymodel.cpp)
target_link_libraries(vc-generate-citymodel pugixml proj shp las)
install(TARGETS vc-generate-citymodel DESTINATION bin)
