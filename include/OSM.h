// OSM I/O
// Anders Logg 2019

#ifndef VC_OSM_H
#define VC_OSM_H

#include <iostream>
#include <pugixml.hpp>

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

    };

};

}

#endif
