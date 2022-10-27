set(PROGRAM run-unittests)

#MESSAGE("Fetching FTP data for unittests")
#execute_process(COMMAND bash -c "${CMAKE_SOURCE_DIR}/data/dtcc-download-data-unittests"
#                OUTPUT_VARIABLE outVar)
#execute_process(COMMAND bash -c "cd ${CMAKE_SOURCE_DIR}/unittests && wget ftps://cloud.dtcc.chalmers.se/Unittests/Unittests.zip -nc --no-check-certificate && unzip -o Unittests.zip && rm -rf Unittests.zip")
#MESSAGE("${outVar}")

add_executable(${PROGRAM} ${PROGRAM}.cpp)
target_link_libraries(${PROGRAM} ${PUGIXML_LIBRARIES} ${PNG_LIBRARIES} ${LIBGEOS_LIBRARIES} shp uuid las dtccproto ${PROTOBUF_LIBRARY} Threads::Threads)
install(TARGETS ${PROGRAM} DESTINATION bin)
