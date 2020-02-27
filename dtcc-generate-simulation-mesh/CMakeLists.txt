set(PROGRAM dtcc-generate-simulation-mesh)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM} dolfin triangle)
install(TARGETS ${PROGRAM} DESTINATION bin)
