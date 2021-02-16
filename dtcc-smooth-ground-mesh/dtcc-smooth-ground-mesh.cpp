// Copyright (C) 2020-2021 Anders Logg, Anton J Olsson
// Licensed under the MIT License

#include <iostream>

#include "CommandLine.h"
#include "JSON.h"
#include "Parameters.h"
#include "Surface.h"
#include "Timer.h"
#include "VTK.h"
#include "VertexSmoother.h"

using namespace DTCC;

void Help()
{
  Error("Usage: dtcc-smooth-surface-mesh Parameters.json");
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
  Point2D origin(parameters.X0, parameters.Y0);
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
  JSON::Write(groundMesh,
              dataDirectory + "VisualizationSmoothedGroundMesh.json", origin);
  VTK::Write(groundMesh, dataDirectory + "VisualizationSmoothedGroundMesh.vtu");

  // Report timings
  Timer::Report("dtcc-smooth-ground-mesh");

  return 0;
}
