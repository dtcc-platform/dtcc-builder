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
#include "Parameters.h"
#include "Polygon.h"

using namespace DTCC;

void Help()
{
  Error("Usage: vc-generate-citymodel Parameters.json");
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
  JSON::Read(cityModel, dataDirectory + "CityModel.json");
  Info(cityModel);

  // Read height map data
  GridField2D heightMap;
  JSON::Read(heightMap, dataDirectory + "HeightMap.json");
  Info(heightMap);

  // Simplify city model
  CityModelGenerator::SimplifyCityModel(cityModel, heightMap,
                                        parameters.MinimalBuildingDistance);
  Info(cityModel);

  // Write city model to file
  JSON::Write(cityModel, dataDirectory + "SimplifiedCityModel.json");

  return 0;
}
