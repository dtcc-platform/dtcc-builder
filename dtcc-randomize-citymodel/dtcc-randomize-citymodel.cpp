// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

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
  Info(parameters);

  // Get parameters
  const std::string dataDirectory = parameters.DataDirectory + "/";

  // Read elevation model
  GridField2D dtm;
  JSON::Read(dtm, dataDirectory + "/DTM.json");
  Info(dtm);

  // Randomize city model
  CityModel cityModel;
  CityModelGenerator::RandomizeCityModel(cityModel, dtm,
                                         parameters.NumRandomBuildings);
  Info(cityModel);
  JSON::Write(cityModel, dataDirectory + "/CityModelRandom.json");

  // Clean city model
  CityModelGenerator::CleanCityModel(cityModel,
                                     parameters.MinimalVertexDistance);
  JSON::Write(cityModel, dataDirectory + "/CityModelClean.json");

  // Simplify city model
  CityModelGenerator::SimplifyCityModel(cityModel,
                                        parameters.MinimalBuildingDistance,
                                        parameters.MinimalVertexDistance);
  JSON::Write(cityModel, dataDirectory + "/CityModelSimple.json");

  return 0;
}
