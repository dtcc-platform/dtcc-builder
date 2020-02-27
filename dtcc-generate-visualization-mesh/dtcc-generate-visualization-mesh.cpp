// vc-generate-visualization-mesh
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

using namespace DTCC;

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

  // Generate 3D surfaces
  std::vector<Surface3D> surfaces = MeshGenerator::GenerateSurfaces3D(
      cityModel, heightMap,
      heightMap.XMin, heightMap.YMin, heightMap.XMax, heightMap.YMax,
      parameters.MeshResolution, parameters.FlatGround);

  // Extract ground surface
  Surface3D groundSurface = surfaces[0];

  // Extract building surfaces
  std::vector<Surface3D> buildingSurfaces;
  for (size_t i = 1; i < surfaces.size(); i++)
    buildingSurfaces.push_back(surfaces[i]);

  // Temporary hack to displace ground mesh
  //for (size_t i = 0; i < groundSurface.Points.size(); i++)
  //  groundSurface.Points[i].z -= 5.0;

  // Convert to FEniCS meshs
  dolfin::Mesh groundMesh;
  dolfin::Mesh buildingMesh;
  FEniCS::ConvertMesh(groundSurface, groundMesh);
  FEniCS::ConvertMesh(buildingSurfaces, buildingMesh);

  // Write to files
  dolfin::File(dataDirectory + "GroundMesh.pvd") << groundMesh;
  dolfin::File(dataDirectory + "BuildingMesh.pvd") << buildingMesh;

  return 0;
}
