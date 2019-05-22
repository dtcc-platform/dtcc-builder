add_executable(vc-generate-heightmap vc-generate-heightmap.cpp)
target_link_libraries(vc-generate-heightmap ${LIBLAS_LIBRARIES})
install(TARGETS vc-generate-heightmap DESTINATION bin)
