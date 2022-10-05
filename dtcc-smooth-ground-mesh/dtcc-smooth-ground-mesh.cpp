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
  error("Usage: dtcc-smooth-surface-mesh Parameters.json");
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
  Parameters p;
  JSON::Read(p, argv[1]);
  info(p);

  // Get data directory
  std::string dataDirectory = p["DataDirectory"];
  dataDirectory += "/";

  // Get origin
  Point2D origin(p["X0"], p["Y0"]);

  // Read ground mesh
  Surface3D groundMesh;
  JSON::Read(groundMesh, dataDirectory + "GroundMesh.json");
  info(groundMesh);

  // Smooth ground mesh
  VertexSmoother::SmoothSurface(groundMesh, p["GroundSmoothing"]);

  // Write to file
  JSON::Write(groundMesh,
              dataDirectory + "VisualizationSmoothedGroundMesh.json", origin);
  VTK::Write(groundMesh, dataDirectory + "VisualizationSmoothedGroundMesh.vtu");

  // Report timings and parameters
  Timer::Report("dtcc-smooth-ground-mesh", dataDirectory);
  info(p);

  return 0;
}
