project(protobuf-bench)
cmake_minimum_required(VERSION 3.16)
set(CMAKE_EXE_LINKER_FLAGS "-static -pthread")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

set(PROGRAM protobuf-bench)
add_executable(${PROGRAM} ${PROGRAM}.cpp dtcc.pb.cc)
target_link_libraries(${PROGRAM} -lprotobuf -p)
install(TARGETS ${PROGRAM} DESTINATION bin)
