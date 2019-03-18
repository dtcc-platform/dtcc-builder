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
        parameters.XMin = json["XMin"];
        parameters.YMin = json["YMin"];
        parameters.XMax = json["XMax"];
        parameters.YMax = json["YMax"];
        parameters.DomainHeight = json["DomainHeight"];
        parameters.HeightMapResolution = json["HeightMapResolution"];
        parameters.MeshResolution = json["MeshResolution"];
    };

    // Write parameters to JSON file
    static void Write(const Parameters& parameters, std::string fileName)
    {
        std::cout << "JSON: " << "Writing parameters to file "
                  << fileName << std::endl;

        // Generate JSON data
        nlohmann::json json;
        json["XMin"] = parameters.XMin;
        json["YMin"] = parameters.YMin;
        json["XMax"] = parameters.XMax;
        json["YMax"] = parameters.YMax;
        json["DomainHeight"] = parameters.DomainHeight;
        json["HeightMapResolution"] = parameters.HeightMapResolution;
        json["MeshResolution"] = parameters.MeshResolution;

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
        heightMap.XMin = json["XMin"];
        heightMap.YMin = json["YMin"];
        heightMap.XMax = json["XMax"];
        heightMap.YMax = json["YMax"];
        heightMap.SizeX = json["SizeX"];
        heightMap.SizeY = json["SizeY"];
        const auto jsonGridData = json["GridData"];
        heightMap.GridData.resize(jsonGridData.size());
        for (size_t i = 0; i < jsonGridData.size(); i++)
            heightMap.GridData[i] = jsonGridData[i];
    };

    // Write height map to JSON file
    static void Write(const HeightMap& heightMap, std::string fileName)
    {
        std::cout << "JSON: " << "Writing height map to file "
                  << fileName << std::endl;

        // Generate JSON data
        nlohmann::json json;
        json["XMin"] = heightMap.XMin;
        json["XMax"] = heightMap.YMin;
        json["YMin"] = heightMap.XMax;
        json["YMax"] = heightMap.YMax;
        json["SizeX"] = heightMap.SizeX;
        json["SizeY"] = heightMap.SizeY;
        json["GridData"] = heightMap.GridData;

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
