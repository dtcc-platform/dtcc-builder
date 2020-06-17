// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#include <iostream>
#include <string>
#include <vector>

#include "CityModel.h"
#include "CityModelGenerator.h"
#include "CommandLine.h"
#include "GridField.h"
#include "JSON.h"
#include "Logging.h"
#include "Parameters.h"
#include "Polygon.h"
#include "SHP.h"

using namespace DTCC;

void Help() { Error("Usage: vc-generate-citymodel Parameters.json"); }

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
  std::vector<Polygon> footprints;
  SHP::Read(footprints, dataDirectory + "PropertyMap.shp");

  // Read height map data
  GridField2D heightMap;
  JSON::Read(heightMap, dataDirectory + "HeightMap.json");
  Info(heightMap);

  // Get origin
  Vector2D origin(parameters.X0, parameters.Y0);

  // Generate city model and transform to new origin
  CityModel cityModel;
  CityModelGenerator::GenerateCityModel(cityModel, footprints, origin,
                                        heightMap.Grid.BoundingBox);

  // Write raw city model to file
  JSON::Write(cityModel, dataDirectory + "CityModelRaw.json");

  // Clean city model and add building heights
  CityModelGenerator::CleanCityModel(cityModel);
  CityModelGenerator::ComputeHeights(cityModel, heightMap);
  Info(cityModel);

  // Write city model to file
  JSON::Write(cityModel, dataDirectory + "CityModel.json");

  // Simplify city model and add building heights
  CityModelGenerator::SimplifyCityModel(cityModel,
                                        parameters.MinimalBuildingDistance);
  CityModelGenerator::ComputeHeights(cityModel, heightMap);
  Info(cityModel);

  // Write simplified city model to file
  JSON::Write(cityModel, dataDirectory + "CityModelSimple.json");

  return 0;
}
