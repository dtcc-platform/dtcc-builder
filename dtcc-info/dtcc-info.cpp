// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#include <iostream>

#include "JSON.h"
#include "LAS.h"
#include "CommandLine.h"
#include "Logging.h"

using namespace DTCC;

void help() { std::cerr << "Usage: vc-info Data.[json,las]" << std::endl; }

int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (argc != 2)
  {
    help();
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
    {
      Parameters parameters;
      JSON::Deserialize(parameters, json);
      Info(parameters);
    }
    else if (typeName == "BoundingBox2D")
    {
      BoundingBox2D boundingBox;
      JSON::Deserialize(boundingBox, json);
      Info(boundingBox);
    }
    else if (typeName == "BoundingBox3D")
    {
      BoundingBox3D boundingBox;
      JSON::Deserialize(boundingBox, json);
      Info(boundingBox);
    }
    else if (typeName == "Grid2D")
    {
      Grid2D grid;
      JSON::Deserialize(grid, json);
      Info(grid);
    }
    else if (typeName == "Grid3D")
    {
      Grid3D grid;
      JSON::Deserialize(grid, json);
      Info(grid);
    }
    else if (typeName == "Mesh2D")
    {
      Mesh2D mesh;
      JSON::Deserialize(mesh, json);
      Info(mesh);
    }
    else if (typeName == "Mesh3D")
    {
      Mesh3D mesh;
      JSON::Deserialize(mesh, json);
      Info(mesh);
    }
    else if (typeName == "GridField2D")
    {
      GridField2D field;
      JSON::Deserialize(field, json);
      Info(field);
    }
    else if (typeName == "GridField3D")
    {
      GridField3D field;
      JSON::Deserialize(field, json);
      Info(field);
    }
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
