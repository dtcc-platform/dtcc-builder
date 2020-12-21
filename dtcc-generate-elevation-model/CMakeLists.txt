set(PROGRAM dtcc-generate-elevation-model)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM} ${LIBLAS_LIBRARIES} ${VTK_LIBRARIES})
install(TARGETS ${PROGRAM} DESTINATION bin)
