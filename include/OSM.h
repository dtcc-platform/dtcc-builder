// OSM I/O
// Anders Logg 2019
// Licensed under the MIT License

#ifndef DTCC_OSM_H
#define DTCC_OSM_H

#include <iostream>
#include <map>
#include <pugixml.hpp>
#include <strings.h>

#include "Vector.h"
#include "Polygon.h"
#include "Logging.h"

namespace DTCC
{

class OSM
{
public:
  // Read polygons from OSM file
  static void Read(std::vector<Polygon> polygons, const std::string& fileName)
  {
    info("OSM: Reading polygons from file " + fileName);

    // Read XML data from file
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(fileName.c_str());

    // Check results
    if (!result)
      throw std::runtime_error(result.description());

    // Extract OSM data
    pugi::xml_node osm = doc.child("osm");

    // Create map for ref --> node
    std::map<size_t, Vector2D> nodeMap;

    // Read nodes (vertices)
    for (pugi::xml_node_iterator it = osm.begin(); it != osm.end(); ++it)
    {
      // Skip if not a node
      if (strcasecmp(it->name(), "node") != 0)
        continue;

      // Add to map
      const size_t ref = it->attribute("id").as_uint();
      const double x = it->attribute("lon").as_double();
      const double y = it->attribute("lat").as_double();
      nodeMap.insert(std::make_pair(ref, Vector2D(x, y)));
    }

    // Diagnostic printing
    Progress("OSM: Parsed " + str(nodeMap.size()) + " nodes");

    // Read ways (footprints)
    for (pugi::xml_node_iterator it = osm.begin(); it != osm.end(); ++it)
    {
      // Skip if not a way
      if (strcasecmp(it->name(), "way") != 0)
        continue;

      // Get building data
      bool isBuilding = false;
      std::string buildingName = "Unknown";
      for (pugi::xml_node_iterator jt = it->begin(); jt != it->end(); ++jt)
      {
        // Skip if not a tag
        if (strcasecmp(jt->name(), "tag") != 0)
          continue;

        // Check whether it is a building
        if (strcasecmp(jt->attribute("k").as_string(), "building") == 0)
          isBuilding = true;

        // Get name of building
        if (strcasecmp(jt->attribute("k").as_string(), "name") == 0)
          buildingName = jt->attribute("v").as_string();
      }

      // Skip if not a building
      if (!isBuilding)
        continue;

      // Create empty polygon
      Polygon polygon;

      // Get vertices
      for (pugi::xml_node_iterator jt = it->begin(); jt != it->end(); ++jt)
      {
        // Skip if not a node
        if (strcasecmp(jt->name(), "nd") != 0)
          continue;

        // Get node and point
        const size_t ref = jt->attribute("ref").as_uint();
        const auto n = nodeMap.find(ref);
        if (n == nodeMap.end())
          throw std::runtime_error("Missing node reference for building.");
        Vector2D p = n->second;

        // Add to polygon
        polygon.Points.push_back(p); // TODO: No Points in Polygon
      }

      // Add polygon
      polygons.push_back(polygon);
    }
  }
};

} // namespace DTCC

#endif
