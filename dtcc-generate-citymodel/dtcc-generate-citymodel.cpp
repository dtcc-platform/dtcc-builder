// vc-generate-citymodel
// Anders Logg 2019

#include <iostream>
#include <string>
#include <vector>

#include "CityModel.h"
#include "CityModelGenerator.h"
#include "CommandLine.h"
#include "HeightMap.h"
#include "JSON.h"
#include "Parameters.h"
#include "Polygon.h"
#include "SHP.h"

using namespace DTCC;

void Help()
{
  std::cerr << "Usage: vc-generate-citymodel Parameters.json" << std::endl;
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

  // Read property map data
  std::vector<Polygon> polygons;
  SHP::Read(polygons, dataDirectory + "PropertyMap.shp");

  // Read height map data
  HeightMap heightMap;
  JSON::Read(heightMap, dataDirectory + "HeightMap.json");
  std::cout << heightMap << std::endl;

  // Generate city model
  CityModel cityModel;
  CityModelGenerator::GenerateCityModel(
      cityModel, polygons, heightMap, parameters.X0, parameters.Y0,
      parameters.XMin, parameters.YMin, parameters.XMax, parameters.YMax);
  std::cout << cityModel << std::endl;

  // Write city model to file
  JSON::Write(cityModel, dataDirectory + "CityModel.json");

  return 0;
}
