// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

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
  Error("Usage: vc-smooth-surface-mesh Parameters.json");
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
  Info(parameters);

  // Get data directory (add trailing slash just in case)
  const std::string dataDirectory = parameters.DataDirectory + "/";

  // Read ground mesh
  Surface3D groundMesh;
  JSON::Read(groundMesh, dataDirectory + "GroundMesh.json");
  Info(groundMesh);

  // Smooth ground mesh
  VertexSmoother::SmoothSurface(groundMesh, parameters.GroundSmoothing);

  // Write to file
  JSON::Write(groundMesh, dataDirectory + "SmoothedGroundMesh.json");
  VTK::Write(groundMesh, dataDirectory + "SmoothedGroundMesh.vtu");

  return 0;
}
