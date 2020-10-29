// Copyright (C) 2020 Anton J Olsson
// Licensed under the MIT LicenseS

#include "CommandLine.h"
#include "Logging.h"
#include <Polygon.h>
#include <Road.h>
#include <SHP.h>
#include <iostream>
#include <json.hpp>
#include <map>

using namespace DTCC;

void Help()
{
  std::cerr << "Usage: vc-generate-roadnetwork fileName.shp" << std::endl;
}

std::vector<Road> MakeRoadObjects(const std::vector<Polygon> &polygons,
                                  const json &attributes);

int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (argc != 2)
  {
    Help();
    return 1;
  }

  // Get .shp filename
  std::string shpFilename = std::string(argv[1]);
  Info(shpFilename);

  // Read road network data
  std::vector<Polygon> roadNetwork;
  nlohmann::json attributes;
  SHP::Read(roadNetwork, shpFilename, &attributes);

  if (attributes.size() != roadNetwork.size())
    throw std::runtime_error("Different number of roads and attribute sets.");
  std::vector<Road> roads = MakeRoadObjects(roadNetwork, attributes);

  return 0;
}
std::vector<Road> MakeRoadObjects(const std::vector<Polygon> &polygons,
                                  const json &attributes)
{
  std::vector<Road> roads;
  for (size_t i = 0; i < polygons.size(); ++i)
  {
    Road road;
    road.Code = attributes[i]["KKOD"];
    road.Category = attributes[i]["KATEGORI"];
    road.Vertices = polygons[i].Vertices;
    roads.push_back(road);
  }
  return roads;
}
