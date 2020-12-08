// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#include <dolfin.h>

#include "CityModel.h"
#include "CommandLine.h"
#include "FEniCS.h"
#include "GridField.h"
#include "JSON.h"
#include "LaplacianSmoother.h"
#include "Logging.h"
#include "Mesh.h"
#include "MeshGenerator.h"
#include "MeshProcessor.h"
#include "Parameters.h"
#include "VTK.h"
#include "VertexSmoother.h"

using namespace DTCC;

void Help()
{
  Error("Usage: dtcc-generate-simulation-mesh Parameters.json");
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

  // Get parameters
  const std::string dataDirectory = parameters.DataDirectory + "/";
  const double h{parameters.MeshResolution};
  const double H{parameters.DomainHeight};

  // Read elevation model (only DTM is used)
  GridField2D dtm;
  JSON::Read(dtm, dataDirectory + "/DTM.json");
  Info(dtm);

  // Smooth elevation model
  VertexSmoother::SmoothField(dtm, parameters.GroundSmoothing);

  // Read city model
  CityModel cityModel;
  JSON::Read(cityModel, dataDirectory + "CityModelSimple.json");
  Info(cityModel);

  // Step 3.1: Generate 2D mesh
  Mesh2D mesh2D;
  MeshGenerator::GenerateMesh2D(mesh2D, cityModel, dtm.Grid.BoundingBox, h);
  Info(mesh2D);
  if (parameters.Debug)
  {
    VTK::Write(mesh2D, dataDirectory + "Step31Mesh.vtu");
  }

  // Step 3.2: Generate 3D mesh (full)
  Mesh3D mesh3D;
  const size_t numLayers = MeshGenerator::GenerateMesh3D(mesh3D, mesh2D, H, h);
  Info(mesh3D);
  if (parameters.Debug)
  {
    Surface3D surf3D;
    MeshProcessor::ExtractBoundary3D(surf3D, mesh3D);
    VTK::Write(mesh3D, dataDirectory + "Step32Mesh.vtu");
    VTK::Write(surf3D, dataDirectory + "Step32Surf.vtu");
  };

  // Step 3.3: Smooth 3D mesh (apply DTM to groundf)
  const double topHeight = dtm.Mean() + H;
  LaplacianSmoother::SmoothMesh3D(mesh3D, cityModel, dtm, topHeight, false);
  Info(mesh3D);
  if (parameters.Debug)
  {
    Surface3D surf3D;
    MeshProcessor::ExtractBoundary3D(surf3D, mesh3D);
    VTK::Write(mesh3D, dataDirectory + "Step33Mesh.vtu");
    VTK::Write(surf3D, dataDirectory + "Step33Surf.vtu");
  }

  // Step 3.4: Trim 3D mesh (remove tets inside buildings)
  MeshGenerator::TrimMesh3D(mesh3D, mesh2D, cityModel, numLayers);
  Info(mesh3D);
  if (parameters.Debug)
  {
    Surface3D surf3D;
    MeshProcessor::ExtractBoundary3D(surf3D, mesh3D);
    VTK::Write(mesh3D, dataDirectory + "Step34Mesh.vtu");
    VTK::Write(surf3D, dataDirectory + "Step34Surf.vtu");
  }

  // Step 3.5: Smooth 3D mesh (apply DTM to ground)
  LaplacianSmoother::SmoothMesh3D(mesh3D, cityModel, dtm, topHeight, true);
  Info(mesh3D);
  if (parameters.Debug)
  {
    Surface3D surf3D;
    MeshProcessor::ExtractBoundary3D(surf3D, mesh3D);
    VTK::Write(mesh3D, dataDirectory + "Step35Mesh.vtu");
    VTK::Write(surf3D, dataDirectory + "Step35Surf.vtu");
  }

  // Write final 3D mesh
  JSON::Write(mesh3D, dataDirectory + "Mesh3D.json");

  // Report timings
  Timer::Report("dtcc-generate-simulation-mesh");

  return 0;
}
