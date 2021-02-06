// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#include <string>
#include <vector>

#include "CommandLine.h"
#include "GridField.h"
#include "JSON.h"
#include "Logging.h"
#include "Parameters.h"
#include "Polygon.h"
#include "SHP.h"
#include "Timer.h"
#include "datamodel//CityModelGenerator.h"
#include "datamodel/CityModel.h"

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

  // Get parameters
  const std::string dataDirectory = parameters.DataDirectory + "/";
  const Vector2D p0{parameters.X0, parameters.Y0};

  // Read property map data
  std::vector<Polygon> footprints;
  std::vector<std::string> UUIDs;
  std::vector<int> entityIDs;
  SHP::Read(footprints, dataDirectory + "/PropertyMap.shp", &UUIDs, &entityIDs);

  // Read elevation models
  GridField2D dsm;
  GridField2D dtm;
  JSON::Read(dsm, dataDirectory + "/DSM.json");
  JSON::Read(dtm, dataDirectory + "/DTM.json");
  Info(dsm);
  Info(dtm);

  // Generate city model and transform to new origin
  CityModel cityModel;
  CityModelGenerator::GenerateCityModel(cityModel, footprints, UUIDs, entityIDs,
                                        p0, dsm.Grid.BoundingBox);
  Info(cityModel);
  JSON::Write(cityModel, dataDirectory + "/CityModelRaw.json");

  // Clean city model and compute building heights
  CityModelGenerator::CleanCityModel(cityModel,
                                     parameters.MinimalVertexDistance);
  CityModelGenerator::ComputeBuildingHeights(cityModel, dsm, dtm);
  Info(cityModel);
  JSON::Write(cityModel, dataDirectory + "/CityModelClean.json");

  // Simplify city model and compute building heights
  CityModelGenerator::SimplifyCityModel(cityModel,
                                        parameters.MinimalBuildingDistance,
                                        parameters.MinimalVertexDistance);
  CityModelGenerator::ComputeBuildingHeights(cityModel, dsm, dtm);
  Info(cityModel);
  JSON::Write(cityModel, dataDirectory + "/CityModelSimple.json");

  // Report timings
  Timer::Report("dtcc-generate-citymodel");

  return 0;
}
