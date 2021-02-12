set(PROGRAM dtcc-generate-mesh)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM} dolfin triangle ${VTK_LIBRARIES})
install(TARGETS ${PROGRAM} DESTINATION bin)
