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
#include "Mesh.h"
#include "MeshGenerator.h"
#include "LaplacianSmoother.h"
#include "Parameters.h"

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
  JSON::Read(heightMap, dataDirectory + "HeightMap.json");
  Info(heightMap);

  // Generate 2D mesh
  Mesh2D mesh2D = MeshGenerator::GenerateMesh2D(
      cityModel, heightMap.Grid.BoundingBox, parameters.MeshResolution);
  Info(mesh2D);

  // Compute ground elevation
  const double groundElevation = heightMap.Min();

  // Generate 3D mesh (excluding height map)
  Mesh3D mesh3D = MeshGenerator::GenerateMesh3D(
      mesh2D, cityModel, groundElevation, parameters.DomainHeight,
      parameters.MeshResolution);
  Info(mesh3D);

  // Convert to FEniCS meshes
  dolfin::Mesh _mesh2D, _mesh3D;
  FEniCS::ConvertMesh(mesh2D, _mesh2D);
  FEniCS::ConvertMesh(mesh3D, _mesh3D);

  // Apply mesh smoothing to account for height map
  LaplacianSmoother::SmoothMesh(_mesh3D, heightMap, cityModel, mesh3D.DomainMarkers,
                                parameters.MeshResolution);

  // Generate height map function (used only for testing/visualization)
  auto z = LaplacianSmoother::GenerateHeightMapFunction(_mesh2D, heightMap);

  // Generate mesh boundary (used only for testing/visualization)
  dolfin::BoundaryMesh _boundary3D(_mesh3D, "exterior");

  // Write mesh to files¨
  std::cout << "vc-generate-mesh: Writing to files..." << std::endl;
  dolfin::File(dataDirectory + "Mesh2D.pvd") << _mesh2D;
  dolfin::File(dataDirectory + "Mesh3D.pvd") << _mesh3D;
  dolfin::File(dataDirectory + "MeshBoundary.pvd") << _boundary3D;
  dolfin::File(dataDirectory + "HeightMap.pvd") << *z;

  // Report timings
  Timer::Report("dtcc-generate-simulation-mesh");

  return 0;
}
