set(BINDINGS _dtcc_builder)
set(pybind11_DIR /usr/lib/cmake/pybind11)
find_package(pybind11 REQUIRED)

pybind11_add_module(${BINDINGS} ${CMAKE_CURRENT_SOURCE_DIR}/cpp/dtcc_builder.cpp)
target_link_libraries(${BINDINGS} PRIVATE triangle ${LIBGEOS_LIBRARIES})

install(TARGETS ${BINDINGS} DESTINATION dtcc_builder)
