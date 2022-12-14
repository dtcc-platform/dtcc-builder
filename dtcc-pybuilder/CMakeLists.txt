set(PYBUILDER _pybuilder)
set(PROGRAM dtcc-build)

find_package(Python COMPONENTS Interpreter REQUIRED)
set(pybind11_DIR /usr/lib/cmake/pybind11)
find_package(pybind11 REQUIRED)
set(DTCC_PROTOBUF_DIR ${CMAKE_SOURCE_DIR}/dtcc-model/protobuf)
set (DTCC_PROTOBUF_FILE ${DTCC_PROTOBUF_DIR}/dtcc.proto)

link_libraries(dtccproto ${PROTOBUF_LIBRARY} Threads::Threads)
pybind11_add_module(${PYBUILDER} ${CMAKE_SOURCE_DIR}/dtcc-pybuilder/dtcc-pybindings/src/pybind_builder.cpp)
target_link_libraries(${PYBUILDER} PRIVATE dolfin triangle ${VTK_LIBRARIES} ${LIBGEOS_LIBRARIES} ${ASSIMP_LIBRARIES} uuid stdc++fs)
set_target_properties(${PYBUILDER} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/dtcc-pybuilder/dtcc/builder)
install(CODE "execute_process(COMMAND ${Python_EXECUTABLE} -m pip install ${CMAKE_CURRENT_LIST_DIR})")
install(FILES ${CMAKE_CURRENT_LIST_DIR}/${PROGRAM} 
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
        GROUP_READ GROUP_EXECUTE 
        WORLD_READ WORLD_EXECUTE 
        DESTINATION bin)