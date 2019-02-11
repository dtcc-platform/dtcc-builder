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

        // Read data from file
        std::ifstream f(fileName);
        nlohmann::json json;
        f >> json;

        // Extract JSON data
        const auto GridData = json["GridData"];
        const auto Width = json["Width"];
        const auto Height = json["Height"];
        const auto A = json["GridMap"]["A"];
        const auto B = json["GridMap"]["B"];
        const auto D = json["GridMap"]["D"];
        const auto E = json["GridMap"]["E"];
        const auto C = json["GridMap"]["C"];
        const auto F = json["GridMap"]["F"];

        // Set height map data
        heightMap.GridData.resize(GridData.size());
        for (size_t i = 0; i < GridData.size(); i++)
            heightMap.GridData[i] = GridData[i];
        heightMap.Width = Width;
        heightMap.Height = Height;
        heightMap.GridMap.A = A;
        heightMap.GridMap.B = B;
        heightMap.GridMap.D = D;
        heightMap.GridMap.E = E;
        heightMap.GridMap.C = C;
        heightMap.GridMap.F = F;
    };

    // Write height map to JSON file
    static void Write(const HeightMap& heightMap, std::string fileName)
    {
        std::cout << "JSON: " << "Writing height map to file "
                  << fileName << std::endl;

        // Generate JSON data
        nlohmann::json json;
        json["Width"] = heightMap.Width;
        json["Height"] = heightMap.Height;
        json["GridData"] = heightMap.GridData;
        json["GridMap"]["A"] = heightMap.GridMap.A;
        json["GridMap"]["B"] = heightMap.GridMap.B;
        json["GridMap"]["D"] = heightMap.GridMap.D;
        json["GridMap"]["E"] = heightMap.GridMap.E;
        json["GridMap"]["C"] = heightMap.GridMap.C;
        json["GridMap"]["F"] = heightMap.GridMap.F;

        // Write to file
        std::ofstream f(fileName);
        f << json;
    }

};

}

#endif
