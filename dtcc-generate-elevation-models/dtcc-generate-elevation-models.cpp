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

void Help() { Error("Usage: dtcc-generate-elevation-models Parameters.json"); }

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

  // Read point cloud data
  PointCloud pointCloud;
  for (auto const &f : CommandLine::ListDirectory(dataDirectory))
  {
    if (CommandLine::EndsWith(f, ".las") || CommandLine::EndsWith(f, ".laz"))
    {
      LAS::Read(pointCloud, dataDirectory + f);
      Info(pointCloud);
    }
  }

  // Set domain size
  double xMin{}, yMin{}, xMax{}, yMax{};
  if (parameters.AutoDomain)
  {
    Progress("Automatically determining domain size:");
    xMin = pointCloud.BoundingBox.P.x - parameters.X0;
    yMin = pointCloud.BoundingBox.P.y - parameters.Y0;
    xMax = pointCloud.BoundingBox.Q.x - parameters.X0;
    yMax = pointCloud.BoundingBox.Q.y - parameters.Y0;
    Progress("  XMin: " + str(pointCloud.BoundingBox.P.x) + " --> " +
             str(xMin));
    Progress("  YMin: " + str(pointCloud.BoundingBox.P.y) + " --> " +
             str(yMin));
    Progress("  XMax: " + str(pointCloud.BoundingBox.Q.x) + " --> " +
             str(xMax));
    Progress("  YMax: " + str(pointCloud.BoundingBox.Q.y) + " --> " +
             str(yMax));
  }
  else
  {
    xMin = parameters.XMin;
    yMin = parameters.YMin;
    xMax = parameters.XMax;
    yMax = parameters.YMax;
  }

  // Generate DSM (including buildings and other objects)
  GridField2D dsm{};
  ElevationModelGenerator::GenerateElevationModel(
      dsm, pointCloud, parameters.X0, parameters.Y0, xMin, yMin, xMax, yMax,
      parameters.ElevationModelResolution);
  Info(dsm);
  JSON::Write(dsm, dataDirectory + "DSM.json");

  // FIXME: We should be able to filter the point cloud (not read again).

  // Read in only ground points
  pointCloud.clear();
  for (auto const &f : CommandLine::ListDirectory(dataDirectory))
  {
    if (CommandLine::EndsWith(f, ".las") || CommandLine::EndsWith(f, ".laz"))
    {
      LAS::Read(pointCloud, dataDirectory + f,
                {2, 9}); // only ground and water points
      Info(pointCloud);
    }
  }

  // Generate DTM (excluding buildings and other objects)
  GridField2D dtm{};
  ElevationModelGenerator::GenerateElevationModel(
      dtm, pointCloud, parameters.X0, parameters.Y0, xMin, yMin, xMax, yMax,
      parameters.ElevationModelResolution);
  Info(dtm);
  JSON::Write(dtm, dataDirectory + "DTM.json");

  // Report timings
  Timer::Report("dtcc-generate-elevation-models");

  return 0;
}
