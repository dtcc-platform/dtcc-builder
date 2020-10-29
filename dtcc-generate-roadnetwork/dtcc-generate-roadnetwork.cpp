// Copyright (C) 2020 Anton J Olsson
// Licensed under the MIT LicenseS

#include "CommandLine.h"
#include "Logging.h"
#include <JSON.h>
#include <Polygon.h>
#include <RoadNetwork.h>
#include <SHP.h>
#include <iostream>
#include <json.hpp>
#include <map>

using namespace DTCC;

void Help()
{
  std::cerr << "Usage: vc-generate-roadnetwork fileName.shp" << std::endl;
}

std::vector<RoadNetwork> GetRoadObjects(const std::vector<Polygon> &polygons,
                                        const json &attributes);

json GetJSONObjects(const std::vector<RoadNetwork> &roads);
std::string GetDataDirectory(const std::string &filename);

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
    throw std::runtime_error("Differing number of roads and attribute sets.");
  std::vector<RoadNetwork> roads = GetRoadObjects(roadNetwork, attributes);

  std::string dataDirectory = GetDataDirectory(shpFilename);
  json jsonNetwork = GetJSONObjects(roads);
  Info(jsonNetwork.dump(4));
  // JSON::Write(jsonNetwork, "RoadNetwork.json");

  return 0;
}
std::string GetDataDirectory(const std::string &filename)
{
  int endPos = filename.find_last_of('/');
  return filename.substr(0, endPos + 1);
}

json GetJSONObjects(const std::vector<RoadNetwork> &roads)
{
  json jsonNetwork = json::array();
  for (const auto &road : roads)
  {
    json jsonRoad = json({});
    JSON::Serialize(road, jsonRoad);
    jsonNetwork.push_back(jsonRoad);
  }
  return jsonNetwork;
}

std::vector<RoadNetwork> GetRoadObjects(const std::vector<Polygon> &polygons,
                                        const json &attributes)
{
  std::vector<RoadNetwork> roads;
  for (size_t i = 0; i < polygons.size(); ++i)
  {
    RoadNetwork road;
    road.Code = attributes[i]["KKOD"];
    road.Category = attributes[i]["KATEGORI"];
    road.Vertices = polygons[i].Vertices;
    roads.push_back(road);
  }
  return roads;
}
