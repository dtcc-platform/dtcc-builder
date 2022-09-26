#include "dtcc.pb.h"
#include "Mesh.h"
#include "Timer.h"

#include <iostream>

void DTCCCreate(DTCC::Mesh3D& mesh)
{
  DTCC::Timer timer("DTCC create");
}

void DTCCAccess(const DTCC::Mesh3D& mesh)
{
  DTCC::Timer timer("DTCC access");
}

void PROTOCreate(PROTO::Mesh3D& mesh)
{
  DTCC::Timer timer("PROTO create");
}

void PROTOAccess(const PROTO::Mesh3D& mesh)
{
  DTCC::Timer timer("PROTO access");
}

int main()
{
  // Run benchmarks for DTCC data structures
  {
    DTCC::Mesh3D mesh{};
    DTCCCreate(mesh);
    DTCCAccess(mesh);
  }

  // Run benchmarks for Protobuf data structures
  {
    PROTO::Mesh3D mesh{};
    PROTOCreate(mesh);
    PROTOAccess(mesh);
  }

  // Report timings
  DTCC::Timer::Report("protobuf-bench");

  return 0;
}
