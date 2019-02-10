// JSON I/O
// Anders Logg 2019

#ifndef JSON_H
#define JSON_H

#include <fstream>
#include <iostream>
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
        std::cout << "JSON: " << "Reading height map from file "
                  << fileName << std::endl;

        std::ifstream f(fileName);
        nlohmann::json json;
        f >> json;

        std::cout << "Not implemented" << std::endl;
    };

    // Write height map to JSON file
    static void Write(const HeightMap& heightMap, std::string fileName)
    {
        std::cout << "JSON: " << "Writing height map to file "
                  << fileName << std::endl;

        std::cout << "Not implemented" << std::endl;
    }

};

}

#endif
