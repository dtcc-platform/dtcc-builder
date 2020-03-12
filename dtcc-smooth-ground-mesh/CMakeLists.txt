set(PROGRAM dtcc-smooth-ground-mesh)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM})
install(TARGETS ${PROGRAM} DESTINATION bin)
