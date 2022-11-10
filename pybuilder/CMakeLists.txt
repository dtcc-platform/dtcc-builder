set(PYBUILDER _pybuilder)
set(PROGRAM dtcc-generate-citymodel.py)

find_package(Python COMPONENTS Interpreter REQUIRED)

set(pybind11_DIR /usr/lib/cmake/pybind11)
find_package(pybind11 REQUIRED)

pybind11_add_module(${PYBUILDER} ${CMAKE_SOURCE_DIR}/pybuilder/pybindings/src/pybind_builder.cpp)
target_link_libraries(${PYBUILDER} PRIVATE dolfin triangle ${VTK_LIBRARIES} ${LIBGEOS_LIBRARIES} ${ASSIMP_LIBRARIES} shp las uuid stdc++fs dtccproto ${PROTOBUF_LIBRARY} Threads::Threads)
set_target_properties(${PYBUILDER} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/pybuilder/pybuilder)
install(CODE "execute_process(COMMAND ${Python_EXECUTABLE} -m pip install ${CMAKE_CURRENT_LIST_DIR})")
install(FILES ${CMAKE_CURRENT_LIST_DIR}/${PROGRAM} DESTINATION bin)