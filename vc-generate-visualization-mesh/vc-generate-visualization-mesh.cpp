// VirtualCity@Chalmers: vc-generate-visualization-mesh
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
  std::cerr << "Usage: vc-generate-visualization-mesh Parameters.json"
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
  JSON::Read(cityModel, dataDirectory + "CityModel.json");
  std::cout << cityModel << std::endl;

  // Read height map data
  HeightMap heightMap;
  JSON::Read(heightMap, dataDirectory + "HeightMap.json");
  std::cout << heightMap << std::endl;

  std::cout << "CHECK: 000" << std::endl;

  // Generate 3D surfaces
  std::vector<Surface3D> surfaces = MeshGenerator::GenerateSurfaces3D(
      cityModel, heightMap, parameters.XMin, parameters.YMin, parameters.XMax,
      parameters.YMax, parameters.MeshResolution, parameters.FlatGround);

  // Convert to FEniCS surface mesh
  dolfin::Mesh surfaceMesh;
  FEniCS::ConvertMesh(surfaces, surfaceMesh);

  // Write to files
  dolfin::File(dataDirectory + "SurfaceMesh.pvd") << surfaceMesh;

  return 0;
}
