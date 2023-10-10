set(BINDINGS _dtcc_builder)

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cpp/)

pybind11_add_module(${BINDINGS} ${SRC_DIR}/dtcc_builder.cpp)

# Find OpenMP
find_package(OpenMP)

if(OpenMP_FOUND)
    message(STATUS "OpenMP found")
    target_compile_definitions(${BINDINGS} PUBLIC _OPENMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    # OpenMP is not available
    message(WARNING "OpenMP not found, disabling OpenMP support")
endif()

target_compile_definitions(${BINDINGS} PRIVATE -DTRILIBRARY -DANSI_DECLARATORS)

if(WIN32)
    target_compile_definitions(${BINDINGS} PRIVATE -DNO_TIMER)
endif()

install(TARGETS ${BINDINGS} DESTINATION dtcc_builder)
