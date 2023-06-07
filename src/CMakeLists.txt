set(BINDINGS _dtcc_builder)
set(pybind11_DIR /usr/lib/cmake/pybind11)
find_package(pybind11 REQUIRED)

pybind11_add_module(${BINDINGS} ${CMAKE_CURRENT_SOURCE_DIR}/cpp/dtcc_builder.cpp)
message(STATUS "Triangle include dir: ${Triangle_INCLUDE_DIR}")
message(STATUS "Triangle library: ${Triangle_LIBRARY}")
target_link_libraries(${BINDINGS} PRIVATE ${Triangle_LIBRARY} ${LIBGEOS_LIBRARIES})

install(TARGETS ${BINDINGS} DESTINATION dtcc_builder)
