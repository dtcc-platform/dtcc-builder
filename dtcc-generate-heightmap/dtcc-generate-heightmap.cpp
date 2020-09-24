// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#include <iostream>
#include <string>
#include <vector>

#include "CommandLine.h"
#include "HeightMapGenerator.h"
#include "JSON.h"
#include "LAS.h"
#include "Logging.h"
#include "Parameters.h"

using namespace DTCC;

void Help()
{
  std::cerr << "Usage: vc-generate-heightmap Parameters.json" << std::endl;
}

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
    Progress("  XMin: " + str(pointCloud.BoundingBox.P.x) + " --> " + str(xMin));
    Progress("  YMin: " + str(pointCloud.BoundingBox.P.y) + " --> " + str(yMin));
    Progress("  XMax: " + str(pointCloud.BoundingBox.Q.x) + " --> " + str(xMax));
    Progress("  YMax: " + str(pointCloud.BoundingBox.Q.y) + " --> " + str(yMax));
  }
  else
  {
    xMin = parameters.XMin;
    yMin = parameters.YMin;
    xMax = parameters.XMax;
    yMax = parameters.YMax;
  }

  // Generate height map
  GridField2D heightMap;
  HeightMapGenerator::GenerateHeightMap(heightMap, pointCloud,
                                        parameters.X0, parameters.Y0,
                                        xMin, yMin, xMax, yMax,
                                        parameters.HeightMapResolution);
  Info(heightMap);

  // Write height map data
  JSON::Write(heightMap, dataDirectory + "HeightMap.json");

  // Report timings
  Timer::Report("dtcc-generate-heightmap");

  return 0;
}
