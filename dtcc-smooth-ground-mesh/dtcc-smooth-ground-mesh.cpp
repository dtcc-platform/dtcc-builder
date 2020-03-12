// vc-smooth-surface-mesh
// Anders Logg 2020

#include <iostream>

#include "CommandLine.h"
#include "Parameters.h"
#include "Surface.h"
#include "VertexSmoother.h"
#include "JSON.h"
#include "VTK.h"

using namespace DTCC;

void Help()
{
  std::cerr << "Usage: vc-smooth-surface-mesh Parameters.json"
            << std::endl;
}

int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (argc != 2)
  {
    Help();
    return 1;
  }

  // Read parameters
  Parameters parameters;
  JSON::Read(parameters, argv[1]);
  std::cout << parameters << std::endl;

  // Get data directory (add trailing slash just in case)
  const std::string dataDirectory = parameters.DataDirectory + "/";

  // Read ground mesh
  Surface3D groundMesh;
  JSON::Read(groundMesh, dataDirectory + "GroundMesh.json");
  std::cout << groundMesh << std::endl;

  // Smooth ground mesh
  VertexSmoother::SmoothSurface(groundMesh, parameters.GroundSmoothing);

  // Write to file
  JSON::Write(groundMesh, dataDirectory + "SmoothedGroundMesh.json");
  VTK::Write(groundMesh, dataDirectory + "SmoothedGroundMesh.vtu");

  return 0;
}
