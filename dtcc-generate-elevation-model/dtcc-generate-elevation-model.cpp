// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#include <iostream>
#include <string>
#include <vector>

#include "CommandLine.h"
#include "ElevationModelGenerator.h"
#include "JSON.h"
#include "LAS.h"
#include "Logging.h"
#include "Parameters.h"
#include "PointCloudProcessor.h"
#include "VTK.h"

using namespace DTCC;

void Help() { Error("Usage: dtcc-generate-elevation-model Parameters.json"); }

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
  const double h{parameters.ElevationModelResolution};
  const Point2D origin{parameters.X0, parameters.Y0};
  Point2D p{parameters.XMin, parameters.YMin};
  Point2D q{parameters.XMax, parameters.YMax};

  // Set bounding box
  p += Vector2D(origin);
  q += Vector2D(origin);
  const BoundingBox2D bbox{p, q};
  Info("Bounding box: " + str(bbox));

  // Read point cloud
  PointCloud pointCloud;
  LAS::ReadDirectory(pointCloud, dataDirectory, bbox);
  pointCloud.SetOrigin(origin);
  Info(pointCloud);

  // Generate DSM (including buildings and other objects)
  GridField2D dsm;
  ElevationModelGenerator::GenerateElevationModel(dsm, pointCloud, h, "DSM");
  Info(dsm);
  JSON::Write(dsm, dataDirectory + "DSM.json");
  if (parameters.Debug)
    VTK::Write(dsm, dataDirectory + "DSM.vts");

  // Filter only ground and water points (color 2 and 9)
  pointCloud = PointCloudProcessor::ClassificationFilter(pointCloud, {2, 9});
  Info(pointCloud);

  // Generate DTM (excluding buildings and other objects)
  GridField2D dtm;
  ElevationModelGenerator::GenerateElevationModel(dtm, pointCloud, h, "DTM");
  Info(dtm);
  JSON::Write(dtm, dataDirectory + "DTM.json");
  if (parameters.Debug)
    VTK::Write(dtm, dataDirectory + "DTM.vts");

  // Report timings
  Timer::Report("dtcc-generate-elevation-models");

  return 0;
}
