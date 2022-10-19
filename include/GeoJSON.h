// Copyright (C) 2020-2021 Anton J Olsson
// Licensed under the MIT License

#ifndef CORE_GEOJSON_H
#define CORE_GEOJSON_H

#include "JSON.h"
#include "Logging.h"
#include "Polyfix.h"
#include "RoadNetworkGenerator.h"
#include <json.hpp>
#include <regex>

using namespace nlohmann;

namespace DTCC_BUILDER
{

/// Class for parsing and converting from GeoJSON format.
class GeoJSON
{

public:
  /// Read GeoJSON file containing RoadNetwork data and write to
  /// RoadNetwork.json. \param filename GeoJSON filename
  static void WriteRoadNetwork(const std::string &filename)
  {
    RoadNetwork roadNetwork = GetRoadNetwork(filename);

    // Compute offset/origin and subtract accordingly
    Point2D origin = BoundingBox2D(roadNetwork.Vertices).P;
    Polyfix::Transform(roadNetwork.Vertices, origin);

    // Change suffix from .geojson to .json
    std::regex reEnding("\\.geo");
    JSON::Write(roadNetwork, std::regex_replace(filename, reEnding, "."),
                origin);
  }

  /// Get RoadNetwork object from GeoJSON file.
  /// \param filename GeoJSON filename
  /// \return RoadNetwork object
  static RoadNetwork GetRoadNetwork(const std::string &filename)
  {
    json geoRoadNetwork, jsonRoadNetwork;
    JSON::Read(geoRoadNetwork, filename);

    return RoadNetworkGenerator::GetRoadNetwork(geoRoadNetwork);
  }
};
} // namespace DTCC_BUILDER

#endif // CORE_GEOJSON_H
