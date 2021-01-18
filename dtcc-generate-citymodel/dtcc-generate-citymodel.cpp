// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#include <string>
#include <vector>

#include "CityModel.h"
#include "CityModelGenerator.h"
#include "CommandLine.h"
#include "GridField.h"
#include "JSON.h"
#include "LAS.h"
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

  // Get parameters
  const std::string dataDirectory{parameters.DataDirectory + "/"};
  const Point2D origin{parameters.X0, parameters.Y0};
  const double minVertexDistance{parameters.MinVertexDistance};

  // Read property map
  std::vector<Polygon> footprints;
  std::vector<std::string> UUIDs;
  std::vector<int> entityIDs;
  SHP::Read(footprints, dataDirectory + "PropertyMap.shp", &UUIDs, &entityIDs);

  // Read point cloud
  PointCloud pointCloud;
  LAS::ReadDirectory(pointCloud, dataDirectory);
  pointCloud.SetOrigin(origin);

  // Read elevation model
  GridField2D dtm;
  JSON::Read(dtm, dataDirectory + "/DTM.json");

  // Generate city model (raw)
  CityModel cityModel;
  CityModelGenerator::GenerateCityModel(cityModel, footprints, UUIDs,
                                        entityIDs);
  cityModel.SetOrigin(origin);
  JSON::Write(cityModel, dataDirectory + "CityModelRaw.json");

  // Clean city model and set heights
  CityModelGenerator::CleanCityModel(cityModel, minVertexDistance);
  CityModelGenerator::ExtractRoofPoints(cityModel, pointCloud);
  CityModelGenerator::ComputeBuildingHeights(cityModel, dtm);
  JSON::Write(cityModel, dataDirectory + "CityModelClean.json");

  // Report timings
  Timer::Report("dtcc-generate-citymodel");

  return 0;
}
