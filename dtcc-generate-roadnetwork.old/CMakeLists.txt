set(PROGRAM dtcc-generate-roadnetwork)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM} shp uuid)
install(TARGETS ${PROGRAM} DESTINATION bin)