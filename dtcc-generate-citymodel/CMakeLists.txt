set(PROGRAM dtcc-generate-citymodel)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM} ${VTK_LIBRARIES} pugixml proj shp las uuid stdc++fs)
install(TARGETS ${PROGRAM} DESTINATION bin)
