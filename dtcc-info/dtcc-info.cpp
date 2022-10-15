// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#include <iostream>

#include "CommandLine.h"
#include "JSON.h"
#include "LAS.h"
#include "Logging.h"
#include "Utils.h"

using namespace DTCCBUILDER;

void Help() { error("Usage: dtcc-info Data.[json,las]"); }

template <class T> void info(nlohmann::json json)
{
  T t;
  JSON::Deserialize(t, json);
  info(t);
}

int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (argc != 2)
  {
    Help();
    return 1;
  }

  // Get filename
  const std::string fileName(argv[1]);

  // Check file type
  if (Utils::EndsWith(fileName, "json"))
  {
    // Read JSON and get type
    nlohmann::json json;
    JSON::Read(json, fileName);
    const std::string typeName = JSON::ToString("Type", json);

    // Check type
    if (typeName == "Parameters")
      info<Parameters>(json);
    else if (typeName == "BoundingBox2D")
      info<BoundingBox2D>(json);
    else if (typeName == "BoundingBox3D")
      info<BoundingBox3D>(json);
    else if (typeName == "Grid2D")
      info<Grid2D>(json);
    else if (typeName == "Grid3D")
      info<Grid3D>(json);
    else if (typeName == "Mesh2D")
      info<Mesh2D>(json);
    else if (typeName == "Mesh3D")
      info<Mesh3D>(json);
    else if (typeName == "GridField2D")
      info<GridField2D>(json);
    else if (typeName == "GridField3D")
      info<GridField3D>(json);
    else if (typeName == "GridVectorField2D")
      info<GridVectorField2D>(json);
    else if (typeName == "GridVectorField3D")
      info<GridVectorField3D>(json);
    else if (typeName == "CityModel")
      info<CityModel>(json);
    else if (typeName == "RoadNetwork")
      info<RoadNetwork>(json);
    else
    {
      error("Unknown JSON type: '" + typeName + "'");
    }
  }
  else if (Utils::EndsWith(fileName, "las"))
  {
    PointCloud pointCloud;
    LAS::Read(pointCloud, fileName);
    info(pointCloud);
  }
  else
  {
    error("Unhandled file type.");
  }

  return 0;
}
