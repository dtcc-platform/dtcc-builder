set(PROGRAM dtcc-generate-cityjson-from-mesh)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM} ${LIBLAS_LIBRARIES} ${ASSIMP_LIBRARIES})
install(TARGETS ${PROGRAM} DESTINATION bin)
