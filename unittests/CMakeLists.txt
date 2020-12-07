set(PROGRAM run-unittests)
add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM} ${PUGIXML_LIBRARIES} ${PNG_LIBRARIES} shp)
install(TARGETS ${PROGRAM} DESTINATION bin)
