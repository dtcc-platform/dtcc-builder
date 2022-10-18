// Copyright (C) 2020 Anton J Olsson
// Licensed under the MIT License

#ifndef CORE_ROADNETWORKGENERATOR_H
#define CORE_ROADNETWORKGENERATOR_H

#include "Hashing.h"
#include "Polygon.h"
#include "SHP.h"
#include "datamodel/RoadNetwork.h"
#include <nlohmann/json.hpp>

using namespace nlohmann;

namespace DTCC_BUILDER
{
/// Class for generating RoadNetwork objects, given .json or .geojson files.
class RoadNetworkGenerator
{
  /// Add edge values (doubles) and/or properties (strings) to RoadNetwork.
  /// \param network RoadNetwork object
  /// \param attributes Values and properties fetched from file
  /// \param attributeIndex Index of edge to add attributes to
  static void AddEdgeProperties(RoadNetwork &network,
                                const json &attributes,
                                size_t attributeIndex)
  {
    json jsonAttributes = attributes["attributes"][attributeIndex];
    for (json::iterator elem = jsonAttributes.begin();
         elem != jsonAttributes.end(); ++elem)
    {
      std::string key = elem.key();
      if (elem.value().is_number_float())
      {
        network.EdgeValues[elem.key()].push_back(elem.value());
        continue;
      }
      std::string value = elem.value().is_string()
                              ? elem.value().get<std::string>()
                              : elem.value().dump();
      network.EdgeProperties[elem.key()].push_back(value);
    }
  }

  /// Get RoadNetwork object given polygons and attributes.
  /// \param polygons Polygon vector describing the road network
  /// \param attributes Attributes for each road edge
  /// \return RoadNetwork object
  static RoadNetwork GetRoadNetwork(const std::vector<Polygon> &polygons,
                                    const json &attributes)
  {
    if (!attributes.is_null() &&
        attributes["polyToAttribute"].size() != polygons.size())
      throw std::runtime_error("Differing number of roads and attribute sets.");

    // Create road edges by fetching unique vertices and pairing their indices
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

          if (!attributes.is_null())
            AddEdgeProperties(network, attributes,
                              attributes["polyToAttribute"][i]);
        }
        vertexIndex++;
      }
    }
    return network;
  }

public:
  /// Get RoadNetwork object given GeoJSON data as nlohmann::json object.
  /// \param geoJson GeoJSON nlohmann::json object
  /// \return RoadNetwork object
  static RoadNetwork GetRoadNetwork(const json &geoJson)
  {
    std::vector<Polygon> polygons;
    json attributes = nullptr;
    for (size_t i = 0; i < geoJson["features"].size(); ++i)
    {
      auto feature = geoJson["features"][i];
      Polygon edge;
      json jsonEdge = feature["geometry"]["coordinates"];
      Point2D p1(jsonEdge[0][0].get<double>(), jsonEdge[0][1].get<double>());
      Point2D p2(jsonEdge[1][0].get<double>(), jsonEdge[1][1].get<double>());
      edge.Vertices = {p1, p2};
      polygons.push_back(edge);

      for (json::iterator property = feature["properties"].begin();
           property != feature["properties"].end(); ++property)
        attributes["attributes"][i][property.key()] = property.value();
      attributes["polyToAttribute"][i] = i;
    }
    return GetRoadNetwork(polygons, attributes);
  }

  /// Get RoadNetwork object from RoadNetwork.json.
  /// \param jsonFilename Name of RoadNetwork JSON file
  /// \return RoadNetwork object
  static RoadNetwork GetRoadNetwork(const std::string &jsonFilename)
  {
    std::vector<Polygon> polygons;
    json attributes;
    SHP::Read(polygons, jsonFilename, nullptr, nullptr, &attributes);
    return GetRoadNetwork(polygons, attributes);
  }
};

} // namespace DTCC_BUILDER

#endif // CORE_ROADNETWORKGENERATOR_H
