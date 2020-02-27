set(PROGRAM vc-info)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM} las)
install(TARGETS ${PROGRAM} DESTINATION bin)
