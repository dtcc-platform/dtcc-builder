set(PYBUILDER _pybuilder)
set(pybind11_DIR /usr/lib/cmake/pybind11)
find_package(pybind11 REQUIRED)


pybind11_add_module(${PYBUILDER} ${CMAKE_CURRENT_SOURCE_DIR}/cpp/builder_bindings.cpp)
target_link_libraries(${PYBUILDER} PRIVATE dolfin triangle ${LIBGEOS_LIBRARIES} stdc++fs)

install(TARGETS ${PYBUILDER} DESTINATION dtcc_builder)