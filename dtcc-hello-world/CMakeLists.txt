set(PROGRAM dtcc-hello-world)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM} ${LIBLAS_LIBRARIES})
install(TARGETS ${PROGRAM} DESTINATION bin)
