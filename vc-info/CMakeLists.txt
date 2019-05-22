add_executable(vc-info vc-info.cpp)
target_link_libraries(vc-info las)
install(TARGETS vc-info DESTINATION bin)
