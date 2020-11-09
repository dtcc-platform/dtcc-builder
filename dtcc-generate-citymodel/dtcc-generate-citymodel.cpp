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
#include "Timer.h"

using namespace DTCC;

void Help() { Error("Usage: dtcc-generate-citymodel Parameters.json"); }

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

  // Read property map data
  std::vector<Polygon> footprints;
  SHP::Read(footprints, parameters.DataDirectory + "/PropertyMap.shp");

  // Read height map data
  GridField2D heightMap;
  JSON::Read(heightMap, parameters.DataDirectory + "/HeightMap.json");
  Info(heightMap);

  // Generate city model and transform to new origin
  Vector2D origin(parameters.X0, parameters.Y0);
  CityModel cityModel;
  CityModelGenerator::GenerateCityModel(cityModel, footprints, origin,
                                        heightMap.Grid.BoundingBox);

  // Write raw city model to file
  JSON::Write(cityModel, parameters.DataDirectory + "/CityModelRaw.json");

  // Clean city model and add building heights
  CityModelGenerator::CleanCityModel(cityModel,
                                     parameters.MinimalVertexDistance);
  CityModelGenerator::ComputeBuildingHeights(cityModel, heightMap);
  Info(cityModel);

  // Write city model to file
  JSON::Write(cityModel, parameters.DataDirectory + "/CityModel.json");

  // Simplify city model and add building heights
  CityModelGenerator::SimplifyCityModel(cityModel,
                                        parameters.MinimalBuildingDistance,
                                        parameters.MinimalVertexDistance);
  CityModelGenerator::ComputeBuildingHeights(cityModel, heightMap);
  Info(cityModel);

  // Write simplified city model to file
  JSON::Write(cityModel, parameters.DataDirectory + "/CityModelSimple.json");

  // Report timings
  Timer::Report("dtcc-generate-citymodel");

  return 0;
}
