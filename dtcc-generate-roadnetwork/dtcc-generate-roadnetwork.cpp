// Copyright (C) 2020 Anton J Olsson
// Licensed under the MIT LicenseS

#include "CommandLine.h"
#include "Logging.h"
#include <Hashing.h>
#include <JSON.h>
#include <Polygon.h>
#include <SHP.h>
#include <datamodel/RoadNetwork.h>
#include <iostream>
#include <json.hpp>
#include <map>

using namespace DTCC;

void Help()
{
  std::cerr << "Usage: dtcc-generate-roadnetwork fileName.shp" << std::endl;
}

RoadNetwork GetRoadNetwork(const std::vector<Polygon> &polygons,
                           json &attributes);
std::string GetDataDirectory(const std::string &filename);
void AddEdgeProperties(RoadNetwork &network,
                       json &attributes,
                       size_t attributeIndex);

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
  std::vector<Polygon> polygons;
  nlohmann::json attributes;
  SHP::Read(polygons, shpFilename, nullptr, nullptr, &attributes);

  if (attributes["polyToAttribute"].size() != polygons.size())
    throw std::runtime_error("Differing number of roads and attribute sets.");
  RoadNetwork network = GetRoadNetwork(polygons, attributes);

  // Serialize and write JSON file
  std::string dataDirectory = GetDataDirectory(shpFilename);
  json jsonNetwork;
  JSON::Serialize(network, jsonNetwork);
  JSON::Write(jsonNetwork, dataDirectory + "RoadNetwork.json");

  return 0;
}

std::string GetDataDirectory(const std::string &filename)
{
  int endPos = filename.find_last_of('/');
  return filename.substr(0, endPos + 1);
}

// Construct unique vertices and edges with indices to the vertices
RoadNetwork GetRoadNetwork(const std::vector<Polygon> &polygons,
                           json &attributes)
{
  RoadNetwork network;
  std::unordered_map<size_t, size_t> hashToVertexIndex;
  size_t vertexIndex = 0;
  for (size_t i = 0; i < polygons.size(); i++)
  {
    size_t numVertices = polygons[i].Vertices.size();
    for (size_t j = 0; j < numVertices; j++)
    {
      size_t hash = Hashing::Hash(polygons[i].Vertices[j]);
      size_t uniqueVertexIndex = network.Vertices.size();
      if (hashToVertexIndex.count(hash) > 0)
        uniqueVertexIndex = hashToVertexIndex[hash];
      else
        network.Vertices.push_back(polygons[i].Vertices[j]);
      hashToVertexIndex.insert(std::make_pair(hash, uniqueVertexIndex));

      if (j > 0)
        network.Edges.back().second = uniqueVertexIndex;
      if (j < numVertices - 1)
      {
        network.Edges.emplace_back(uniqueVertexIndex, -1);
        AddEdgeProperties(network, attributes,
                          attributes["polyToAttribute"][i]);
      }
      vertexIndex++;
    }
  }
  return network;
}

void AddEdgeProperties(RoadNetwork &network,
                       json &attributes,
                       size_t attributeIndex)
{
  json jsonAttributes = attributes["attributes"][attributeIndex];
  for (json::iterator elem = jsonAttributes.begin();
       elem != jsonAttributes.end(); ++elem)
  {
    std::string key = elem.key();
    std::string value = elem.value().is_string()
                            ? elem.value().get<std::string>()
                            : elem.value().dump();
    network.EdgeProperties[elem.key()].push_back(value);
  }
}
