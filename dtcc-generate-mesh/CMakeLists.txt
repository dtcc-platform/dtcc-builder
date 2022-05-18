set(PROGRAM dtcc-generate-mesh)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM} dolfin triangle geos_c ${VTK_LIBRARIES} shp las uuid)
install(TARGETS ${PROGRAM} DESTINATION bin)
