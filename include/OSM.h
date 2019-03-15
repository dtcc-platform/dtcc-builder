// OSM I/O
// Anders Logg 2019

#ifndef VC_OSM_H
#define VC_OSM_H

#include <iostream>
#include <strings.h>
#include <map>
#include <pugixml.hpp>

#include "Point.h"
#include "Building.h"
#include "CoordinateSystem.h"

namespace VirtualCity
{

class OSM
{
public:

    // Read city model from OSM file
    static void Read(CityModel& cityModel, std::string fileName)
    {
        std::cout << "OSM: " << "Reading city model from file "
                  << fileName << std::endl;

        // Read XML data from file
        pugi::xml_document doc;
        pugi::xml_parse_result result = doc.load_file(fileName.c_str());

        // Check results
        if (!result)
            throw std::runtime_error(result.description());

        // Extract OSM data
        pugi::xml_node osm = doc.child("osm");

        // Create map for ref --> node
        std::map<size_t, Point2D> nodeMap;

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
            nodeMap.insert(std::make_pair(ref, Point2D(x, y)));
        }

        // Diagnostic printing
        std::cout << "OSM: Parsed " << nodeMap.size() << " nodes" << std::endl;

        // Read ways (footprints)
        for (pugi::xml_node_iterator it = osm.begin(); it != osm.end(); ++it)
        {
            // Skip if not a way
            if (strcasecmp(it->name(), "way") != 0)
                continue;

            // Get building data
            bool isBuilding = false;
            std::string buildingName = "Unknown";
            for (pugi::xml_node_iterator jt = it->begin();
                    jt != it->end(); ++jt)
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

            // Create empty building
            Building building;

            // Read footprint
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
                Point2D p = n->second;

                // Add to building footprint
                building.Footprint.push_back(p);
            }

            // Add building to city model
            cityModel.Buildings.push_back(building);

            // Diagnostic printing
            std::cout << "OSM: Added building " << buildingName << " ("
                      << building.Footprint.size() << " nodes)" << std::endl;
        }
    }

};

}

#endif
