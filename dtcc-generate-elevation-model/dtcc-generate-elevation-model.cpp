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

  // Get data directory (add trailing slash just in case)
  const std::string dataDirectory = parameters.DataDirectory + "/";

  // Read point cloud data (all *.las and *.laz files in data directory)
  PointCloud pointCloud;
  for (auto const &f : CommandLine::ListDirectory(dataDirectory))
  {
    if (CommandLine::EndsWith(f, ".las") || CommandLine::EndsWith(f, ".laz"))
    {
      LAS::Read(pointCloud, dataDirectory + f);
      Info(pointCloud);
    }
  }

  // Get domain size, either from point cloud bounding box or parameters
  BoundingBox2D bbox;
  if (parameters.AutoDomain)
  {
    Progress("Automatically determining domain size:");
    bbox.P.x = pointCloud.BoundingBox.P.x - parameters.X0;
    bbox.P.y = pointCloud.BoundingBox.P.y - parameters.Y0;
    bbox.Q.x = pointCloud.BoundingBox.Q.x - parameters.X0;
    bbox.Q.y = pointCloud.BoundingBox.Q.y - parameters.Y0;
  }
  else
  {
    bbox.P.x = parameters.XMin;
    bbox.P.y = parameters.YMin;
    bbox.Q.x = parameters.XMax;
    bbox.Q.y = parameters.YMax;
  }
  Progress("Domain bounding box: " + str(bbox));

  // Get origin and grid resolution
  const Point2D x0{parameters.X0, parameters.Y0};
  const double h{parameters.ElevationModelResolution};

  // Generate DSM (including buildings and other objects)
  GridField2D dsm;
  ElevationModelGenerator::GenerateElevationModel(dsm, pointCloud, x0, bbox, h);
  Info(dsm);
  JSON::Write(dsm, dataDirectory + "DSM.json");

  // FIXME: We should be able to filter the point cloud (not read again).

  // Read in only ground and water points (color 2 and 9)
  pointCloud.clear();
  for (auto const &f : CommandLine::ListDirectory(dataDirectory))
  {
    if (CommandLine::EndsWith(f, ".las") || CommandLine::EndsWith(f, ".laz"))
    {
      LAS::Read(pointCloud, dataDirectory + f, {2, 9});
      Info(pointCloud);
    }
  }

  // Generate DTM (excluding buildings and other objects)
  GridField2D dtm;
  ElevationModelGenerator::GenerateElevationModel(dtm, pointCloud, x0, bbox, h);
  Info(dtm);
  JSON::Write(dtm, dataDirectory + "DTM.json");

  // Report timings
  Timer::Report("dtcc-generate-elevation-models");

  return 0;
}
