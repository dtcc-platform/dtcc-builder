set(BINDINGS _dtcc_builder)
set(pybind11_DIR /usr/lib/cmake/pybind11)
find_package(pybind11 REQUIRED)

pybind11_add_module(${BINDINGS} ${CMAKE_CURRENT_SOURCE_DIR}/cpp/dtcc_builder.cpp)

target_compile_definitions(${BINDINGS} PRIVATE -DTRILIBRARY -DANSI_DECLARATORS)
if(WIN32)
    target_compile_definitions(${BINDINGS} PRIVATE -DNO_TIMER)
endif()

target_link_libraries(${BINDINGS} PRIVATE ${LIBGEOS_LIBRARIES})

install(TARGETS ${BINDINGS} DESTINATION dtcc_builder)
