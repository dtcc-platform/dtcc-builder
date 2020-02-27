set(PROGRAM vc-simplify-citymodel)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM} pugixml)
install(TARGETS ${PROGRAM} DESTINATION bin)
