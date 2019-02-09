// JSON I/O
// Anders Logg 2019

#ifndef JSON_H
#define JSON_H

#include <fstream>
#include <json.hpp>
#include "HeightMap.h"

namespace VirtualCity
{

class JSON
{
public:

    // Read height map from JSON file
    static void Read(HeightMap& heightMap, std::string fileName)
    {
        std::ifstream f(fileName);
        nlohmann::json json;
        f >> json;
    };

    // Write height map to JSOON file
    static void Write(const HeightMap& heightMap, std::string fileName)
    {
        std::cout << "Not implemented" << std::endl;
    }

};

}

#endif
