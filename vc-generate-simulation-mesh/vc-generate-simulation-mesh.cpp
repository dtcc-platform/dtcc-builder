// VirtualCity@Chalmers: vc-generate-simulation-mesh
// Anders Logg 2018

#include <dolfin.h>
#include <iostream>
#include <string>
#include <vector>

#include "CityModel.h"
#include "CommandLine.h"
#include "FEniCS.h"
#include "HeightMap.h"
#include "JSON.h"
#include "Mesh.h"
#include "MeshGenerator.h"
#include "MeshSmoother.h"
#include "Parameters.h"

using namespace VirtualCity;

void Help()
{
  std::cerr << "Usage: vc-generate-simulation-mesh Parameters.json"
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

  // Read city model data
  CityModel cityModel;
  JSON::Read(cityModel, dataDirectory + "SimplifiedCityModel.json");
  std::cout << cityModel << std::endl;

  // Read height map data
  HeightMap heightMap;
  JSON::Read(heightMap, dataDirectory + "HeightMap.json");
  std::cout << heightMap << std::endl;

  // Generate 2D mesh
  Mesh2D mesh2D = MeshGenerator::GenerateMesh2D(
      cityModel, parameters.XMin, parameters.YMin, parameters.XMax,
      parameters.YMax, parameters.MeshResolution);
  std::cout << mesh2D << std::endl;

  // Compute ground elevation
  const double groundElevation = heightMap.Min();

  // Generate 3D mesh (excluding height map)
  Mesh3D mesh3D = MeshGenerator::GenerateMesh3D(
      mesh2D, cityModel, groundElevation, parameters.DomainHeight,
      parameters.MeshResolution);
  std::cout << mesh3D << std::endl;

  // Convert to FEniCS meshes
  dolfin::Mesh _mesh2D, _mesh3D;
  FEniCS::ConvertMesh(mesh2D, _mesh2D);
  FEniCS::ConvertMesh(mesh3D, _mesh3D);

  // Apply mesh smoothing to account for height map
  MeshSmoother::SmoothMesh(_mesh3D, heightMap, cityModel, mesh3D.DomainMarkers,
                           parameters.MeshResolution);

  // Generate height map function (used only for testing/visualization)
  auto z = MeshSmoother::GenerateHeightMapFunction(_mesh2D, heightMap);

  // Generate mesh boundary (used only for testing/visualization)
  dolfin::BoundaryMesh _boundary3D(_mesh3D, "exterior");

  // Write mesh to files¨
  std::cout << "vc-generate-mesh: Writing to files..." << std::endl;
  dolfin::File(dataDirectory + "Mesh2D.pvd") << _mesh2D;
  dolfin::File(dataDirectory + "Mesh3D.pvd") << _mesh3D;
  dolfin::File(dataDirectory + "MeshBoundary.pvd") << _boundary3D;
  dolfin::File(dataDirectory + "HeightMap.pvd") << *z;

  return 0;
}
