// Copyright (C) 2020-2022 Anders Logg, Anton J Olsson, Dag Wästberg
// Licensed under the MIT License

#include <experimental/filesystem>
#include <string>
#include <vector>

// Needs to come before JSON (nlohmann) include because of sloppy
// namespacing in VTK (typedef detail)...
#include "VTK.h"

// DTCC includes
#include "CityModelGenerator.h"
#include "CommandLine.h"
#include "ElevationModelGenerator.h"
#include "GridField.h"
#include "JSON.h"
#include "LAS.h"
#include "Logging.h"
#include "ParameterProcessor.h"
#include "Parameters.h"
#include "Polygon.h"
#include "SHP.h"
#include "Timer.h"
#include "Utils.h"
#include "VertexSmoother.h"

using namespace DTCC;

void Help() { Error("Usage: dtcc-generate-citymodel Parameters.json"); }

int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (CommandLine::HasOption("-h", argc, argv))
  {
    Help();
    return 0;
  }
  Parameters p = ParameterProcessor::ProcessArgs(argc, argv);

  // Get directories
  const std::string dataDirectory = p["DataDirectory"];
  const std::string outputDirectory = p["OutputDirectory"];
  Info("Loding data from directory: " + dataDirectory);
  Info("Saving data to directory:   " + outputDirectory);

  // Start timer
  Timer timer("Step 1: Generate city model");

  // Read property map
  std::vector<Polygon> footprints;
  std::vector<std::string> UUIDs;
  std::vector<int> entityIDs;

  SHP::Read(footprints, dataDirectory + "PropertyMap.shp", &UUIDs, &entityIDs);
  Info("Loaded " + str(footprints.size()) + " building footprints");

  // Set bounding box
  BoundingBox2D bbox;
  Point2D O;
  if (p["AutoDomain"])
  {
    bbox = BoundingBox2D(footprints, p["DomainMargin"]);
    Info("Bounding box of footprints: " + str(bbox));
    BoundingBox2D lasBBox;
    LAS::BoundsDirectory(lasBBox, dataDirectory);
    Info("Bounding box of point cloud: " + str(lasBBox));
    bbox.Intersect(lasBBox);
    O = bbox.P;
    p["X0"] = O.x;
    p["Y0"] = O.y;
    p["XMin"] = 0.0;
    p["YMin"] = 0.0;
    p["XMax"] = bbox.Q.x - bbox.P.x;
    p["YMax"] = bbox.Q.y - bbox.P.y;
  }
  else
  {
    O = Point2D(p["X0"], p["Y0"]);
    const double xMin = p["XMin"];
    const double xMax = p["XMax"];
    const double yMin = p["YMin"];
    const double yMax = p["YMax"];
    const Point2D P{O.x + xMin, O.y + yMin};
    const Point2D Q{O.x + xMax, O.y + yMax};
    bbox = BoundingBox2D(P, Q);
  }

  // Check size of bounding box
  Info("Bounding box of city model: " + str(bbox));
  if (bbox.Area() < 100.0)
  {
    Error("Domain too small to generate a city model");
    return 1;
  }

  // Read point cloud (only points inside bounding box)
  PointCloud pointCloud;
  LAS::ReadDirectory(pointCloud, dataDirectory, bbox,
                     p["NaiveVegitationFilter"]);

  // Check point cloud
  if (pointCloud.Empty())
    Error("Point cloud is empty. Check LiDaR quality or the X{0,Min,Max}, "
          "Y{0,Min,Max} values in Parameters.json");
  pointCloud.SetOrigin(O);
  pointCloud.BuildHasClassifications();
  Info(pointCloud);

  // Remove outliers from point cloud
  PointCloudProcessor::RemoveOutliers(pointCloud, p["OutlierMargin"]);

  if (p["NaiveVegitationFilter"])
  {
    // Remove vegetation from point cloud
    PointCloudProcessor::NaiveVegetationFilter(pointCloud);
  }

  // Generate DTM (excluding buildings and other objects)
  GridField2D dtm;
  ElevationModelGenerator::GenerateElevationModel(
      dtm, pointCloud, {2, 9}, p["ElevationModelResolution"]);
  Info(dtm);

  // Generate DSM (including buildings and other objects)
  GridField2D dsm;
  ElevationModelGenerator::GenerateElevationModel(
      dsm, pointCloud, {}, p["ElevationModelResolution"]);
  Info(dsm);

  // Smooth elevation model (only done for DTM)
  VertexSmoother::SmoothField(dtm, p["GroundSmoothing"]);

  // Generate raw city model
  CityModel cityModel;

  std::string modelName = p["ModelName"];
  if (modelName.size() == 0)
    modelName = Utils::GetFilename(dataDirectory);

  cityModel.Name = modelName;
  CityModelGenerator::GenerateCityModel(cityModel, footprints, UUIDs, entityIDs,
                                        bbox, p["MinBuildingDistance"],
                                        p["MinBuildingSize"]);
  cityModel.SetOrigin(O);
  Info(cityModel);

  // Clean city model and compute heights
  CityModelGenerator::CleanCityModel(cityModel, p["MinVertexDistance"]);
  CityModelGenerator::ExtractBuildingPoints(
      cityModel, pointCloud, p["GroundMargin"], p["OutlierMargin"]);

  // Remove outliers from roofs using RANSAC
  // Best for very noisy data like old landmäteriet data
  // if buildings are classified the RANSAC outlier remover
  // will almost certainly do more harm than good
  bool classifiedBuildings = pointCloud.HasClassification(6);

  if (!classifiedBuildings && p["RANSACOutlierRemover"])
  {
    CityModelGenerator::BuildingPointsRANSACOutlierRemover(
        cityModel, p["RANSACOutlierMargin"], p["RANSACIterations"]);
  }

  // Remove roof outliers using Statistics Outlier Removal
  if (p["StatisticalOutlierRemover"])
  {
    CityModelGenerator::BuildingPointsOutlierRemover(
        cityModel, p["OutlierNeighbors"], p["OutlierSTD"]);
  }

  CityModelGenerator::ComputeBuildingHeights(
      cityModel, dtm, p["GroundPercentile"], p["RoofPercentile"]);

  // Stop timer
  timer.Stop();

  // Write JSON
  if (p["WriteJSON"])
  {
    JSON::Write(dtm, outputDirectory + "DTM.json", O);
    JSON::Write(dsm, outputDirectory + "DSM.json", O);
    JSON::Write(cityModel, outputDirectory + "CityModel.json", O);
  }

  // Write VTK
  if (p["WriteVTK"])
  {
    VTK::Write(dtm, outputDirectory + "DTM.vts");
    VTK::Write(dsm, outputDirectory + "DSM.vts");
  }

  // Report timings and parameters
  const std::string prefix = outputDirectory + "/dtcc-generate-citymodel";
  Timer::Report("Timings for dtcc-generate-citymodel",
                prefix + "-timings.json");
  JSON::Write(p, prefix + "-parameters.json", 4);
  Info(p);

  return 0;
}
