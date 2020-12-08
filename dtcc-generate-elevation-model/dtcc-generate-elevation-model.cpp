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
  const Point2D p{parameters.XMin, parameters.YMin};
  const Point2D q{parameters.XMax, parameters.YMax};
  const Point2D p0{parameters.X0, parameters.Y0};
  BoundingBox2D bbox{p, q};
  const double h{parameters.ElevationModelResolution};

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

  // Automatically determine domain size if auto
  if (parameters.AutoDomain)
  {
    Progress("Automatically determining domain size");
    bbox = pointCloud.BoundingBox;
    bbox.P -= Vector2D{p0};
    bbox.Q -= Vector2D{p0};
  }
  Progress("Domain bounding box: " + str(bbox));

  // Generate DSM (including buildings and other objects)
  GridField2D dsm;
  ElevationModelGenerator::GenerateElevationModel(dsm, pointCloud, p0, bbox, h);
  Info(dsm);
  JSON::Write(dsm, dataDirectory + "DSM.json");
  // FIXME: Not implemented
  // VTK::Write(dsm, dataDirectory + "DSM.vtu");

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
  ElevationModelGenerator::GenerateElevationModel(dtm, pointCloud, p0, bbox, h);
  Info(dtm);
  JSON::Write(dtm, dataDirectory + "DTM.json");
  // FIXME: Not implemented
  // VTK::Write(dtm, dataDirectory + "DTM.vtu");

  // Report timings
  Timer::Report("dtcc-generate-elevation-models");

  return 0;
}
