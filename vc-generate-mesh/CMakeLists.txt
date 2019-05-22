add_executable(vc-generate-mesh vc-generate-mesh.cpp)
target_link_libraries(vc-generate-mesh dolfin triangle)
install(TARGETS vc-generate-mesh DESTINATION bin)
