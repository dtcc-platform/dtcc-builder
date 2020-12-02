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
  Vector2D origin(parameters.X0, parameters.Y0);

  // Read property map data
  std::vector<Polygon> footprints;
  std::vector<std::string> UUIDs;
  std::vector<int> entityIDs;
  SHP::Read(footprints, UUIDs, entityIDs,
            parameters.DataDirectory + "/PropertyMap.shp");
  // Read elevation models
  GridField2D dsm{};
  GridField2D dtm{};
  JSON::Read(dsm, parameters.DataDirectory + "/DSM.json");
  JSON::Read(dtm, parameters.DataDirectory + "/DTM.json");
  Info(dsm);
  Info(dtm);

  // Generate city model and transform to new origin
  CityModel cityModel{};
  CityModelGenerator::GenerateCityModel(cityModel, footprints, UUIDs, entityIDs,
                                        origin, dsm.Grid.BoundingBox);
  Info(cityModel);
  JSON::Write(cityModel, parameters.DataDirectory + "/CityModelRaw.json");

  // Clean city model and add building heights
  CityModelGenerator::CleanCityModel(cityModel,
                                     parameters.MinimalVertexDistance);
  CityModelGenerator::ComputeBuildingHeights(cityModel, dsm, dtm);
  Info(cityModel);
  JSON::Write(cityModel, parameters.DataDirectory + "/CityModelClean.json");

  // Simplify city model and add building heights
  CityModelGenerator::SimplifyCityModel(cityModel,
                                        parameters.MinimalBuildingDistance,
                                        parameters.MinimalVertexDistance);
  CityModelGenerator::ComputeBuildingHeights(cityModel, dsm, dtm);
  Info(cityModel);
  JSON::Write(cityModel, parameters.DataDirectory + "/CityModelSimple.json");

  // Report timings
  Timer::Report("dtcc-generate-citymodel");

  return 0;
}
