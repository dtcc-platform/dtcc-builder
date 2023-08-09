set(BINDINGS _dtcc_builder)

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cpp/)

pybind11_add_module(${BINDINGS} ${SRC_DIR}/dtcc_builder.cpp)
target_compile_definitions(${BINDINGS}  PRIVATE -DTRILIBRARY -DANSI_DECLARATORS)

message(STATUS "Triangle include dir: ${Triangle_INCLUDE_DIR}")
message(STATUS "Triangle library: ${Triangle_LIBRARY}")
target_link_libraries(${BINDINGS} PRIVATE triangle)

install(TARGETS ${BINDINGS} DESTINATION dtcc_builder)
