set(PROGRAM run-tests)

add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM})
install(TARGETS ${PROGRAM} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)
