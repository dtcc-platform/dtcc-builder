// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#include <iostream>
#include <string>
#include <vector>

#include "CityModel.h"
#include "CityModelGenerator.h"
#include "CommandLine.h"
#include "JSON.h"
#include "GridField.h"
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
  Info(parameters);

  // Get data directory (add trailing slash just in case)
  const std::string dataDirectory = parameters.DataDirectory + "/";

  // Read property map data
  std::vector<Polygon> polygons;
  SHP::Read(polygons, dataDirectory + "PropertyMap.shp");

  // Read height map data
  GridField2D heightMap;
  JSON::Read(heightMap, dataDirectory + "HeightMap.json");
  Info(heightMap);

  // Generate city model
  CityModel cityModel;
  CityModelGenerator::GenerateCityModel(cityModel, polygons, heightMap,
                                        parameters.X0, parameters.Y0,
                                        heightMap.Grid.BoundingBox.P.x,
                                        heightMap.Grid.BoundingBox.P.y,
                                        heightMap.Grid.BoundingBox.Q.x,
                                        heightMap.Grid.BoundingBox.Q.y);
  Info(cityModel);

  // Write city model to file
  JSON::Write(cityModel, dataDirectory + "CityModel.json");

  return 0;
}
