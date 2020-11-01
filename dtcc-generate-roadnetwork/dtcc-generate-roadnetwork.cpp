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

RoadNetwork GetRoadNetwork(const std::vector<Polygon> &polygons,
                           json &attributes);

std::string GetDataDirectory(const std::string &filename);
std::string GetHash(const Point2D &vertex);

void AddEdgeProperties(RoadNetwork &network, json &attributes, size_t polyNum);
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
  std::vector<Polygon> vertices;
  nlohmann::json attributes;
  SHP::Read(vertices, shpFilename, &attributes);

  if (attributes.size() != vertices.size())
    throw std::runtime_error("Differing number of roads and attribute sets.");
  RoadNetwork network = GetRoadNetwork(vertices, attributes);

  std::string dataDirectory = GetDataDirectory(shpFilename);
  json jsonNetwork;
  JSON::Serialize(network, jsonNetwork);
  Info(jsonNetwork.dump(4));

  JSON::Write(jsonNetwork,
              "/home/dtcc/core/data/roadNetworkSample/RoadNetwork.json");

  return 0;
}

std::string GetDataDirectory(const std::string &filename)
{
  int endPos = filename.find_last_of('/');
  return filename.substr(0, endPos + 1);
}

RoadNetwork GetRoadNetwork(const std::vector<Polygon> &polygons,
                           json &attributes)
{
  RoadNetwork network;
  std::unordered_map<std::string, size_t> hashToVertexIndex;
  size_t vertexIndex = 0;
  for (size_t i = 0; i < polygons.size(); i++)
  {
    size_t numVertices = polygons[i].Vertices.size();
    for (size_t j = 0; j < polygons[i].Vertices.size(); j++)
    {
      std::string hash = GetHash(polygons[i].Vertices[j]);

      size_t uniqueVertexIndex = network.Vertices.size();
      if (hashToVertexIndex.count(hash) > 0)
        uniqueVertexIndex = hashToVertexIndex[hash];
      else
        network.Vertices.push_back(polygons[i].Vertices[j]);
      hashToVertexIndex.insert(std::make_pair(hash, uniqueVertexIndex));

      if (j > 0)
        network.Edges.back().second = uniqueVertexIndex;
      if (j < numVertices - 1)
        network.Edges.emplace_back(uniqueVertexIndex, -1);
      AddEdgeProperties(network, attributes, i);
      vertexIndex++;
    }
  }
  return network;
}
void AddEdgeProperties(RoadNetwork &network, json &attributes, size_t polyNum)
{
  for (json::iterator elem = attributes[polyNum].begin();
       elem != attributes[polyNum].end(); ++elem)
  {
    std::string value = elem.value().is_string()
                            ? elem.value().get<std::string>()
                            : elem.value().dump();
    network.EdgeProperties[elem.key()].push_back(value);
  }
}

std::string GetHash(const Point2D &vertex)
{
  // return std::hash<std::string>{}(vertex.__str__());
  return vertex.__str__();
}
