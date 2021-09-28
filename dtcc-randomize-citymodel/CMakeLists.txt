set(PROGRAM dtcc-randomize-citymodel)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM} ${VTK_LIBRARIES})
install(TARGETS ${PROGRAM} DESTINATION bin)
