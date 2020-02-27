set(PROGRAM dtcc-generate-heightmap)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM} ${LIBLAS_LIBRARIES})
install(TARGETS ${PROGRAM} DESTINATION bin)
