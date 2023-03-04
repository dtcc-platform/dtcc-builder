set(PYBUILDER _pybuilder)
set(pybind11_DIR /usr/lib/cmake/pybind11)
find_package(pybind11 REQUIRED)

link_libraries(${DTCCMODEL_LIBRARIES} ${PROTOBUF_LIBRARY} Threads::Threads)
pybind11_add_module(${PYBUILDER} ${CMAKE_CURRENT_SOURCE_DIR}/pybind_builder.cpp)
target_link_libraries(${PYBUILDER} PRIVATE dolfin triangle ${VTK_LIBRARIES} ${LIBGEOS_LIBRARIES} uuid stdc++fs)

install(TARGETS ${PYBUILDER} DESTINATION dtcc_builder)