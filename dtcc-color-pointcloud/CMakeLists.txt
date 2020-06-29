set(PROGRAM dtcc-color-pointcloud)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM} ${GDAL_LIBRARIES} ${LIBLAS_LIBRARIES})
install(TARGETS ${PROGRAM} DESTINATION bin)
