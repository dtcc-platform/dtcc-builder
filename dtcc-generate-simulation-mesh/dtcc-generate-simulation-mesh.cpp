// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#include <dolfin.h>
#include <iostream>
#include <string>
#include <vector>

#include "CityModel.h"
#include "CommandLine.h"
#include "FEniCS.h"
#include "GridField.h"
#include "JSON.h"
#include "LaplacianSmoother.h"
#include "Logging.h"
#include "Mesh.h"
#include "MeshGenerator.h"
#include "Parameters.h"
#include "VertexSmoother.h"

using namespace DTCC;

void Help()
{
  Error("Usage: vc-generate-simulation-mesh Parameters.json");
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

  // Read city model data
  CityModel cityModel;
  JSON::Read(cityModel, dataDirectory + "CityModelSimple.json");
  Info(cityModel);

  // Read height map data
  GridField2D heightMap;
  JSON::Read(heightMap, dataDirectory + "GroundMap.json");
  Info(heightMap);

  // Smooth height map
  VertexSmoother::SmoothField(heightMap, parameters.GroundSmoothing);

  // Compute absolute height of top of domain
  const double topHeight = heightMap.Mean() + parameters.DomainHeight;

  // Generate 2D mesh
  Mesh2D mesh2D;
  MeshGenerator::GenerateMesh2D(mesh2D, cityModel, heightMap.Grid.BoundingBox,
                                parameters.MeshResolution);
  Info(mesh2D);

  // Write mesh for debugging
  if (parameters.Debug)
  {
    dolfin::Mesh _mesh2D;
    FEniCS::ConvertMesh(mesh2D, _mesh2D);
    auto z = LaplacianSmoother::GenerateHeightMapFunction(_mesh2D, heightMap);
    FEniCS::Write(_mesh2D, dataDirectory + "Mesh2D.pvd");
    FEniCS::Write(*z, dataDirectory + "HeightMap.pvd");
  }

  // Generate 3D mesh
  Mesh3D mesh3D;
  const size_t numLayers = MeshGenerator::GenerateMesh3D(
      mesh3D, mesh2D, parameters.DomainHeight, parameters.MeshResolution);
  Info(mesh3D);

  // Uncomment to write to file for debugging
  if (parameters.Debug)
  {
    dolfin::Mesh _mesh3D;
    FEniCS::ConvertMesh(mesh3D, _mesh3D);
    dolfin::BoundaryMesh _boundary3D(_mesh3D, "exterior");
    FEniCS::Write(_mesh3D, dataDirectory + "Mesh3DFull.pvd");
    FEniCS::Write(_boundary3D, dataDirectory + "Surface3DFull.pvd");
  }

  // Apply mesh smoothing to ground
  LaplacianSmoother::SmoothMesh(mesh3D, cityModel, heightMap, topHeight, false);
  Info(mesh3D);

  // Uncomment to write to file for debugging
  if (parameters.Debug)
  {
    dolfin::Mesh _mesh3D;
    FEniCS::ConvertMesh(mesh3D, _mesh3D);
    dolfin::BoundaryMesh _boundary3D(_mesh3D, "exterior");
    FEniCS::Write(_mesh3D, dataDirectory + "Mesh3DFullSmoothed.pvd");
    FEniCS::Write(_boundary3D, dataDirectory + "Surface3DFullSmoothed.pvd");
  }

  // Trim 3D mesh (remove tets inside buildings)
  MeshGenerator::TrimMesh3D(mesh3D, mesh2D, cityModel, numLayers);
  Info(mesh3D);

  // Uncomment to write to file for debugging
  if (parameters.Debug)
  {
    dolfin::Mesh _mesh3D;
    FEniCS::ConvertMesh(mesh3D, _mesh3D);
    dolfin::BoundaryMesh _boundary3D(_mesh3D, "exterior");
    FEniCS::Write(_mesh3D, dataDirectory + "Mesh3DTrimmed.pvd");
    FEniCS::Write(_boundary3D, dataDirectory + "Surface3DTrimmed.pvd");
  }

  // Apply mesh smoothing to ground and buildings
  LaplacianSmoother::SmoothMesh(mesh3D, cityModel, heightMap, topHeight, true);

  // Uncomment to write to file for debugging
  if (parameters.Debug)
  {
    dolfin::Mesh _mesh3D;
    FEniCS::ConvertMesh(mesh3D, _mesh3D);
    dolfin::BoundaryMesh _boundary3D(_mesh3D, "exterior");
    FEniCS::Write(_mesh3D, dataDirectory + "Mesh3DTrimmedSmoothed.pvd");
    FEniCS::Write(_boundary3D, dataDirectory + "Surface3DTrimmedSmoothed.pvd");
  }

  // Convert to FEniCS meshes (used for testing/visualization)
  dolfin::Mesh _mesh2D, _mesh3D;
  FEniCS::ConvertMesh(mesh2D, _mesh2D);
  FEniCS::ConvertMesh(mesh3D, _mesh3D);

  // Write to files
  Info("vc-generate-simulation-mesh: Writing to files...");
  JSON::Write(mesh2D, dataDirectory + "Mesh2D.json");
  JSON::Write(mesh3D, dataDirectory + "Mesh3D.json");

  // Report timings
  Timer::Report("dtcc-generate-simulation-mesh");

  return 0;
}
