set(PROGRAM dtcc-cityjson-hello-world)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM})
install(TARGETS ${PROGRAM} DESTINATION bin)
