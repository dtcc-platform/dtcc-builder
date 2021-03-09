// Copyright (C) 2020-2021 Anders Logg, Anton J Olsson
// Licensed under the MIT License

#include <string>
#include <vector>

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

void Help() { Error("Usage: dtcc-randomize-citymodel Parameters.json"); }

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
  const Point2D origin{parameters.X0, parameters.Y0};
  Info(parameters);

  // Get parameters
  const std::string dataDirectory{parameters.DataDirectory + "/"};
  const double minBuildingDistance{parameters.MinBuildingDistance};
  const double minVertexDistance{parameters.MinVertexDistance};
  const size_t numRandomBuildings = parameters.NumRandomBuildings;

  // Read elevation model
  GridField2D dtm;
  JSON::Read(dtm, dataDirectory + "/DTM.json");
  Info(dtm);

  // Randomize city model
  CityModel cityModel;
  CityModelGenerator::RandomizeCityModel(cityModel, dtm, numRandomBuildings);
  Info(cityModel);
  JSON::Write(cityModel, dataDirectory + "/CityModelRandom.json", origin);

  // Clean city model
  CityModelGenerator::CleanCityModel(cityModel, minVertexDistance);
  JSON::Write(cityModel, dataDirectory + "/CityModelClean.json", origin);

  // Simplify city model
  CityModelGenerator::SimplifyCityModel(cityModel, minBuildingDistance,
                                        minVertexDistance);
  JSON::Write(cityModel, dataDirectory + "/CityModelSimple.json", origin);

  return 0;
}
