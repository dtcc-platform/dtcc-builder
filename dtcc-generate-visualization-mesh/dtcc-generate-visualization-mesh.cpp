// vc-generate-visualization-mesh
// Anders Logg 2018

#include <dolfin.h>
#include <iostream>
#include <string>
#include <vector>

#include "CityModel.h"
#include "CommandLine.h"
#include "HeightMap.h"
#include "JSON.h"
#include "VTK.h"
#include "Mesh.h"
#include "MeshGenerator.h"
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

  // FIXME: Consider moving this somewhere else as a utility for merging surfaces
  // Extract building surface as one common surface
  std::cout << "Merging building meshes..." << std::endl;
  Surface3D buildingSurface;
  size_t numPoints = 0;
  size_t numCells = 0;
  for (size_t i = 1; i < surfaces.size(); i++)
  {
    numPoints += surfaces[i].Points.size();
    numCells += surfaces[i].Cells.size();
  }
  buildingSurface.Points.resize(numPoints);
  buildingSurface.Cells.resize(numCells);
  size_t k = 0;
  size_t l = 0;
  for (size_t i = 1; i < surfaces.size(); i++)
  {
    for (size_t j = 0; j < surfaces[i].Cells.size(); j++)
    {
      Simplex2D c = surfaces[i].Cells[j];
      c.v0 += k;
      c.v1 += k;
      c.v2 += k;
      buildingSurface.Cells[l++] = c;
    }
    for (size_t j = 0; j < surfaces[i].Points.size(); j++)
      buildingSurface.Points[k++] = surfaces[i].Points[j];
  }

  // Write to files
  JSON::Write(groundSurface, dataDirectory + "GroundMesh.json");
  VTK::Write(groundSurface, dataDirectory + "GroundMesh.vtu");
  VTK::Write(buildingSurface, dataDirectory + "BuildingMesh.vtu");

  return 0;
}
