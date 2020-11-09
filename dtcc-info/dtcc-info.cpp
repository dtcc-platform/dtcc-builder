// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#include <iostream>

#include "JSON.h"
#include "LAS.h"
#include "CommandLine.h"
#include "Logging.h"

using namespace DTCC;

void Help()
{
  Error("Usage: dtcc-info Data.[json,las]");
}

template <class T> void Info(nlohmann::json json)
{
  T t;
  JSON::Deserialize(t, json);
  Info(t);
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
  if (CommandLine::EndsWith(fileName, "json"))
  {
    // Read JSON and get type
    nlohmann::json json;
    JSON::Read(json, fileName);
    const std::string typeName = JSON::ToString("Type", json);

    // Check type
    if (typeName == "Parameters")
      Info<Parameters>(json);
    else if (typeName == "BoundingBox2D")
      Info<BoundingBox2D>(json);
    else if (typeName == "BoundingBox3D")
      Info<BoundingBox3D>(json);
    else if (typeName == "Grid2D")
      Info<Grid2D>(json);
    else if (typeName == "Grid3D")
      Info<Grid3D>(json);
    else if (typeName == "Mesh2D")
      Info<Mesh2D>(json);
    else if (typeName == "Mesh3D")
      Info<Mesh3D>(json);
    else if (typeName == "GridField2D")
      Info<GridField2D>(json);
    else if (typeName == "GridField3D")
      Info<GridField3D>(json);
    else if (typeName == "GridVectorField2D")
      Info<GridVectorField2D>(json);
    else if (typeName == "GridVectorField3D")
      Info<GridVectorField3D>(json);
    else if (typeName == "CityModel")
      Info<CityModel>(json);
    else if (typeName == "RoadNetwork")
      Info<Road>(json);
    else
    {
      Error("Unknown JSON type: '" + typeName + "'");
    }
  }
  else if (CommandLine::EndsWith(fileName, "las"))
  {
    PointCloud pointCloud;
    LAS::Read(pointCloud, fileName);
    Info(pointCloud);
  }
  else
  {
    Error("Unhandled file type.");
  }

  return 0;
}
