set(BINDINGS _dtcc_builder)

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cpp/)

pybind11_add_module(${BINDINGS} ${SRC_DIR}/dtcc_builder.cpp)

# Find OpenMP
find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
    include_directories(SYSTEM ${OpenMP_CXX_INCLUDE_DIRS})
    include_directories(SYSTEM ${OpenMP_omp_INCLUDE_DIR})

    target_link_libraries(${BINDINGS} PUBLIC OpenMP::OpenMP_CXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
else()
    # OpenMP is not available
    message(WARNING "OpenMP not found, disabling OpenMP support")
endif()

target_compile_definitions(${BINDINGS} PRIVATE -DTRILIBRARY -DANSI_DECLARATORS)

if(WIN32)
    target_compile_definitions(${BINDINGS} PRIVATE -DNO_TIMER)
endif()

install(TARGETS ${BINDINGS} DESTINATION dtcc_builder)
