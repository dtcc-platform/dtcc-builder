// JSON I/O
// Anders Logg 2019

#ifndef VC_JSON_H
#define VC_JSON_H

#include <fstream>
#include <iostream>
#include <json.hpp>

#include "Parameters.h"
#include "HeightMap.h"
#include "CityModel.h"

namespace VirtualCity
{

class JSON
{
public:

    // Read parameters from JSON file
    static void Read(Parameters& parameters, std::string fileName)
    {
        std::cout << "JSON: " << "Reading parameters from file "
                  << fileName << std::endl;

        // Read data from file
        std::ifstream f(fileName);
        nlohmann::json json;
        f >> json;

        // Extract JSON data
        parameters.DomainCenterX = json["DomainCenterX"];
        parameters.DomainCenterY = json["DomainCenterY"];
        parameters.DomainRadius = json["DomainRadius"];
        parameters.DomainHeight = json["DomainHeight"];
        parameters.MeshSize = json["MeshSize"];
        parameters.HeightMapStride = json["HeightMapStride"];
    };

    // Write parameters to JSON file
    static void Write(const Parameters& parameters, std::string fileName)
    {
        std::cout << "JSON: " << "Writing parameters to file "
                  << fileName << std::endl;

        // Generate JSON data
        nlohmann::json json;
        json["DomainCenterX"] = parameters.DomainCenterX;
        json["DomainCenterY"] = parameters.DomainCenterY;
        json["DomainRadius"] = parameters.DomainRadius;
        json["DomainHeight"] = parameters.DomainHeight;
        json["MeshSize"] = parameters.MeshSize;
        json["HeightMapStride"] = parameters.HeightMapStride;

        // Write to file
        std::ofstream f(fileName);
        f << json;
    }

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
        const auto jsonGridData = json["GridData"];
        heightMap.GridData.resize(jsonGridData.size());
        for (size_t i = 0; i < jsonGridData.size(); i++)
            heightMap.GridData[i] = jsonGridData[i];
        heightMap.Width = json["Width"];
        heightMap.Height = json["Height"];
        heightMap.GridMap.A = json["GridMap"]["A"];
        heightMap.GridMap.B = json["GridMap"]["B"];
        heightMap.GridMap.D = json["GridMap"]["D"];
        heightMap.GridMap.E = json["GridMap"]["E"];
        heightMap.GridMap.C = json["GridMap"]["C"];
        heightMap.GridMap.F = json["GridMap"]["F"];
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

    // Read city model from JSON file
    static void Read(CityModel& cityModel, std::string fileName)
    {
        std::cout << "JSON: " << "Reading city model from file "
                  << fileName << std::endl;

        // Read data from file
        std::ifstream f(fileName);
        nlohmann::json json;
        f >> json;

        // Extract JSON data
        auto jsonBuildings = json;
        cityModel.Buildings.resize(jsonBuildings.size());
        for (size_t i = 0; i < jsonBuildings.size(); i++)
        {
            auto jsonBuilding = jsonBuildings[i];
            auto jsonFootprint = jsonBuilding["Footprint"];
            cityModel.Buildings[i].Footprint.resize(jsonFootprint.size());
            for (size_t j = 0; j < jsonFootprint.size(); j++)
            {
                cityModel.Buildings[i].Footprint[j].x = jsonFootprint[j]["x"];
                cityModel.Buildings[i].Footprint[j].y = jsonFootprint[j]["y"];
            }
            cityModel.Buildings[i].Height = jsonBuilding["Height"];
        }
    }

    // Write city model to JSON file
    static void Write(CityModel& cityModel, std::string fileName)
    {
        std::cout << "JSON: " << "Writing city model to file "
                  << fileName << std::endl;

        // Generate JSON data
        auto jsonBuildings = nlohmann::json::array();
        for (auto const building : cityModel.Buildings)
        {
            auto jsonBuilding = nlohmann::json::object();
            jsonBuilding["Footprint"] = nlohmann::json::array();
            for (auto const point : building.Footprint)
            {
                auto jsonPoint = nlohmann::json::object();
                jsonPoint["x"] = point.x;
                jsonPoint["y"] = point.y;
                jsonBuilding["Footprint"].push_back(jsonPoint);
            }
            jsonBuilding["Height"] = building.Height;
            jsonBuildings.push_back(jsonBuilding);
        }
        nlohmann::json json;
        json = jsonBuildings;

        // Write to file
        std::ofstream f(fileName);
        f << json;
    }

};

}

#endif
