// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#include <dolfin.h>

#include "CommandLine.h"
#include "FEniCS.h"
#include "GridField.h"
#include "JSON.h"
#include "LaplacianSmoother.h"
#include "Logging.h"
#include "Mesh.h"
#include "MeshGenerator.h"
#include "MeshProcessor.h"
#include "OBJ.h"
#include "Parameters.h"
#include "VTK.h"
#include "VertexSmoother.h"
#include "citymodel/CityModel.h"

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
    VTK::Write(mesh2D, dataDirectory + "Step31Mesh.vtu");

  // Step 3.2: Generate 3D mesh (full)
  Mesh3D mesh;
  const size_t numLayers = MeshGenerator::GenerateMesh3D(mesh, mesh2D, H, h);
  Info(mesh);
  if (parameters.Debug)
  {
    Surface3D boundary;
    MeshProcessor::ExtractBoundary3D(boundary, mesh);
    VTK::Write(mesh, dataDirectory + "Step32Mesh.vtu");
    VTK::Write(boundary, dataDirectory + "Step32Boundary.vtu");
  };

  // Step 3.3: Smooth 3D mesh (apply DTM to ground)
  const double topHeight = dtm.Mean() + H;
  LaplacianSmoother::SmoothMesh3D(mesh, cityModel, dtm, topHeight, false);
  Info(mesh);
  if (parameters.Debug)
  {
    Surface3D boundary;
    MeshProcessor::ExtractBoundary3D(boundary, mesh);
    VTK::Write(mesh, dataDirectory + "Step33Mesh.vtu");
    VTK::Write(boundary, dataDirectory + "Step33Boundary.vtu");
  }

  // Step 3.4: Trim 3D mesh (remove tets inside buildings)
  MeshGenerator::TrimMesh3D(mesh, mesh2D, cityModel, numLayers);
  Info(mesh);
  if (parameters.Debug)
  {
    Surface3D boundary;
    MeshProcessor::ExtractBoundary3D(boundary, mesh);
    VTK::Write(mesh, dataDirectory + "Step34Mesh.vtu");
    VTK::Write(boundary, dataDirectory + "Step34Boundary.vtu");
  }

  // Step 3.5: Smooth 3D mesh (apply DTM to ground)
  LaplacianSmoother::SmoothMesh3D(mesh, cityModel, dtm, topHeight, true);
  Info(mesh);
  if (parameters.Debug)
  {
    Surface3D boundary;
    MeshProcessor::ExtractBoundary3D(boundary, mesh);
    VTK::Write(mesh, dataDirectory + "Step35Mesh.vtu");
    VTK::Write(boundary, dataDirectory + "Step35Boundary.vtu");
  }

  // Extract boundary of final mesh
  Surface3D boundary;
  MeshProcessor::ExtractBoundary3D(boundary, mesh);

  // Write final mesh and boundary
  JSON::Write(mesh, dataDirectory + "Mesh.json");
  JSON::Write(boundary, dataDirectory + "Boundary.json");

  // Report timings
  Timer::Report("dtcc-generate-simulation-mesh");

  return 0;
}
