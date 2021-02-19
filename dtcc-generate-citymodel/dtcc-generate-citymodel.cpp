// Copyright (C) 2020-2021 Anders Logg, Anton J Olsson
// Licensed under the MIT License

#include <string>
#include <vector>

// Needs to come before JSON (nlohmann) include because of sloppy
// namespacing in VTK (typedef detail)...
#include "VTK.h"

// DTCC includes
#include "CommandLine.h"
#include "ElevationModelGenerator.h"
#include "GridField.h"
#include "JSON.h"
#include "LAS.h"
#include "Logging.h"
#include "Parameters.h"
#include "Polygon.h"
#include "SHP.h"
#include "Timer.h"
#include "VertexSmoother.h"
#include "datamodel/CityModel.h"
#include "datamodel/CityModelGenerator.h"

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
  Parameters p;
  JSON::Read(p, argv[1]);
  Info(p);

  // Set data directory
  const std::string dataDirectory{p.DataDirectory + "/"};

  // Set bounding box
  const Point2D O{p.X0, p.Y0};
  const Point2D P{p.XMin + p.X0, p.YMin + p.Y0};
  const Point2D Q{p.XMax + p.X0, p.YMax + p.Y0};
  const BoundingBox2D bbox{P, Q};
  Info("Bounding box: " + str(bbox));

  // Read point cloud (only points inside bounding box)
  PointCloud pointCloud;
  LAS::ReadDirectory(pointCloud, dataDirectory, bbox);

  // Check point cloud
  if (pointCloud.Empty())
  {
  Error("Point cloud is empty. Check LiDaR quality or the X{0,Min,Max}, Y{0,Min,Max} values in Parameters.json");
  return 1;
  }
  pointCloud.SetOrigin(O);
  Info(pointCloud);

  // Generate DTM (excluding buildings and other objects)
  GridField2D dtm;
  ElevationModelGenerator::GenerateElevationModel(dtm, pointCloud, {2, 9},
                                                  p.ElevationModelResolution);
  Info(dtm);

  // Generate DSM (including buildings and other objects)
  GridField2D dsm;
  ElevationModelGenerator::GenerateElevationModel(dsm, pointCloud, {},
                                                  p.ElevationModelResolution);
  Info(dsm);

  // Smooth elevation model (only done for DTM)
  VertexSmoother::SmoothField(dtm, p.GroundSmoothing);

  // Read property map
  std::vector<Polygon> footprints;
  std::vector<std::string> UUIDs;
  std::vector<int> entityIDs;
  SHP::Read(footprints, dataDirectory + "PropertyMap.shp", &UUIDs, &entityIDs);

  // Generate raw city model
  CityModel cityModel;
  CityModelGenerator::GenerateCityModel(cityModel, footprints, UUIDs, entityIDs,
                                        bbox);
  cityModel.SetOrigin(O);
  Info(cityModel);

  // Clean city model and compute heights
  CityModelGenerator::CleanCityModel(cityModel, p.MinVertexDistance);
  CityModelGenerator::ExtractBuildingPoints(cityModel, pointCloud,
                                            p.GroundMargin);
  CityModelGenerator::ComputeBuildingHeights(cityModel, dtm, p.GroundPercentile,
                                             p.RoofPercentile);

  // Write to file
  JSON::Write(dtm, dataDirectory + "DTM.json", O);
  JSON::Write(dsm, dataDirectory + "DSM.json", O);
  JSON::Write(cityModel, dataDirectory + "CityModel.json", O);
  if (p.Debug)
  {
    VTK::Write(dtm, dataDirectory + "DTM.vts");
    VTK::Write(dsm, dataDirectory + "DSM.vts");
  }

  // Report timings
  Timer::Report("dtcc-generate-citymodel");

  return 0;
}
