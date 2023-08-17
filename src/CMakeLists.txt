set(BINDINGS _dtcc_builder)

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cpp/)

pybind11_add_module(${BINDINGS} ${SRC_DIR}/dtcc_builder.cpp)

target_compile_definitions(${BINDINGS} PRIVATE -DTRILIBRARY -DANSI_DECLARATORS)
if(WIN32)
    target_compile_definitions(${BINDINGS} PRIVATE -DNO_TIMER)
endif()


install(TARGETS ${BINDINGS} DESTINATION dtcc_builder)
