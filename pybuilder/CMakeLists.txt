set(PYBUILDER _pybuilder)

set(pybind11_DIR /usr/lib/cmake/pybind11)
find_package(pybind11 REQUIRED)

pybind11_add_module(${PYBUILDER} ${CMAKE_SOURCE_DIR}/pybuilder/pybindings/src/pybind_builder.cpp)
target_link_libraries(${PYBUILDER} PRIVATE dolfin triangle ${VTK_LIBRARIES} ${LIBGEOS_LIBRARIES} ${ASSIMP_LIBRARIES} shp las uuid stdc++fs)
set_target_properties(${PYBUILDER} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/pybuilder/pybuilder)