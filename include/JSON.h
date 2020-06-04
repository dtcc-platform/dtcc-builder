// JSON I/O
// Anders Logg 2019

#ifndef DTCC_JSON_H
#define DTCC_JSON_H

#include <fstream>
#include <iostream>
#include <json.hpp>

#include "CityModel.h"
#include "HeightMap.h"
#include "Mesh.h"
#include "Surface.h"
#include "Parameters.h"

namespace DTCC
{

class JSON
{
public:
  // Get type of JSON file
  static std::string ReadType(std::string fileName)
  {
    // Read data from file
    std::ifstream f(fileName);
    nlohmann::json json;
    f >> json;

    // Check file type
    if (json.find("Type") == json.end())
      throw std::runtime_error("Not a DTCC JSON file.");

    return json["Type"];
  }

  // Read parameters from JSON file
  static void Read(Parameters &parameters, std::string fileName)
  {
    std::cout << "JSON: Reading parameters from file " << fileName << std::endl;

    // Read data from file
    std::ifstream f(fileName);
    nlohmann::json json;
    f >> json;

    // Check file type
    if (json.find("Type") == json.end() || json["Type"] != "Parameters")
      throw std::runtime_error("Not a DTCC Parameters JSON file.");

    // Extract JSON data
    parameters.DataDirectory = TryReadString("DataDirectory", json);
    parameters.X0 = TryReadDouble("X0", json);
    parameters.Y0 = TryReadDouble("Y0", json);
    parameters.XMin = TryReadDouble("XMin", json);
    parameters.YMin = TryReadDouble("YMin", json);
    parameters.XMax = TryReadDouble("XMax", json);
    parameters.YMax = TryReadDouble("YMax", json);
    parameters.AutoDomain = TryReadBool("AutoDomain", json);
    parameters.DomainHeight = TryReadDouble("DomainHeight", json);
    parameters.HeightMapResolution = TryReadDouble("HeightMapResolution", json);
    parameters.MeshResolution = TryReadDouble("MeshResolution", json);
    parameters.MinimalBuildingDistance =
        TryReadDouble("MinimalBuildingDistance", json);
    parameters.FlatGround = TryReadBool("FlatGround", json);
    parameters.GroundSmoothing = TryReadInt("GroundSmoothing", json);
  };

  // Write parameters to JSON file
  static void Write(const Parameters &parameters, std::string fileName)
  {
    std::cout << "JSON: Writing parameters to file " << fileName << std::endl;

    // Set JSON data
    nlohmann::json json;
    json["Type"] = "Parameters";
    json["DataDirectory"] = parameters.DataDirectory;
    json["X0"] = parameters.X0;
    json["Y0"] = parameters.Y0;
    json["YMin"] = parameters.YMin;
    json["XMin"] = parameters.XMin;
    json["YMin"] = parameters.YMin;
    json["XMax"] = parameters.XMax;
    json["YMax"] = parameters.YMax;
    json["AutoDomain"] = parameters.AutoDomain;
    json["DomainHeight"] = parameters.DomainHeight;
    json["HeightMapResolution"] = parameters.HeightMapResolution;
    json["MeshResolution"] = parameters.MeshResolution;
    json["MinimalBuildingDistance"] = parameters.MinimalBuildingDistance;
    json["FlatGround"] = parameters.FlatGround;
    json["GroundSmoothings"] = parameters.GroundSmoothing;

    // Write to file
    std::ofstream f(fileName);
    f << json;
  }

  // Read height map from JSON file
  static void Read(HeightMap &heightMap, std::string fileName)
  {
    std::cout << "JSON: Reading height map from file " << fileName << std::endl;

    // Read data from file
    std::ifstream f(fileName);
    nlohmann::json json;
    f >> json;

    // Check file type
    if (json.find("Type") == json.end() || json["Type"] != "HeightMap")
      throw std::runtime_error("Not a DTCC HeightMap JSON file.");

    // Extract JSON data
    heightMap.XMin = json["XMin"];
    heightMap.YMin = json["YMin"];
    heightMap.XMax = json["XMax"];
    heightMap.YMax = json["YMax"];
    heightMap.XSize = json["XSize"];
    heightMap.YSize = json["YSize"];
    heightMap.XStep = json["XStep"];
    heightMap.YStep = json["YStep"];
    const auto jsonGridData = json["GridData"];
    heightMap.GridData.resize(jsonGridData.size());
    for (size_t i = 0; i < jsonGridData.size(); i++)
      heightMap.GridData[i] = jsonGridData[i];
  };

  // Write height map to JSON file
  static void Write(const HeightMap &heightMap, std::string fileName)
  {
    std::cout << "JSON: Writing height map to file " << fileName << std::endl;

    // Set JSON data
    nlohmann::json json;
    json["Type"] = "HeightMap";
    json["XMin"] = heightMap.XMin;
    json["YMin"] = heightMap.YMin;
    json["XMax"] = heightMap.XMax;
    json["YMax"] = heightMap.YMax;
    json["XSize"] = heightMap.XSize;
    json["YSize"] = heightMap.YSize;
    json["XStep"] = heightMap.XStep;
    json["YStep"] = heightMap.YStep;
    json["GridData"] = heightMap.GridData;

    // Write to file
    std::ofstream f(fileName);
    f << json;
  }

  // Read city model from JSON file
  static void Read(CityModel &cityModel, std::string fileName)
  {
    std::cout << "JSON: Reading city model from file " << fileName << std::endl;

    // Read data from file
    std::ifstream f(fileName);
    nlohmann::json json;
    f >> json;

    // Check file type
    if (json.find("Type") == json.end() || json["Type"] != "CityModel")
      throw std::runtime_error("Not a DTCC CityModel JSON file.");

    // Extract JSON data
    auto jsonBuildings = json["Buildings"];
    cityModel.Buildings.resize(jsonBuildings.size());
    for (size_t i = 0; i < jsonBuildings.size(); i++)
    {
      auto jsonBuilding = jsonBuildings[i];
      auto jsonFootprint = jsonBuilding["Footprint"];
      cityModel.Buildings[i].Footprint.Vertices.resize(jsonFootprint.size());
      for (size_t j = 0; j < jsonFootprint.size(); j++)
      {
        cityModel.Buildings[i].Footprint.Vertices[j].x = jsonFootprint[j]["x"];
        cityModel.Buildings[i].Footprint.Vertices[j].y = jsonFootprint[j]["y"];
      }
      cityModel.Buildings[i].Height = jsonBuilding["Height"];
    }
  }

  // Write city model to JSON file
  static void Write(CityModel &cityModel, std::string fileName)
  {
    std::cout << "JSON: Writing city model to file " << fileName << std::endl;

    // Generate JSON data
    auto jsonBuildings = nlohmann::json::array();
    for (auto const &building : cityModel.Buildings)
    {
      auto jsonBuilding = nlohmann::json::object();
      jsonBuilding["Footprint"] = nlohmann::json::array();
      for (auto const &point : building.Footprint.Vertices)
      {
        auto jsonPoint = nlohmann::json::object();
        jsonPoint["x"] = point.x;
        jsonPoint["y"] = point.y;
        jsonBuilding["Footprint"].push_back(jsonPoint);
      }
      jsonBuilding["Height"] = building.Height;
      jsonBuildings.push_back(jsonBuilding);
    }

    // Set JSON data
    nlohmann::json json;
    json["Type"] = "CityModel";
    json["Buildings"] = jsonBuildings;

    // Write to file
    std::ofstream f(fileName);
    f << json;
  }

  // FIXME: I/O for mesh domain markers currently not supported

  // Read 2D mesh from JSON file
  static void Read(Mesh2D &mesh, std::string fileName)
  {
    std::cout << "JSON: Reading 2D mesh from file " << fileName << std::endl;

    // Read data from file
    std::ifstream f(fileName);
    nlohmann::json json;
    f >> json;

    // Check file type
    if (json.find("Type") == json.end() || json["Type"] != "Mesh2D")
      throw std::runtime_error("Not a DTCC Mesh2D JSON file.");

    // Extract point data
    const auto jsonPoints = json["Points"];
    mesh.Vertices.resize(jsonPoints.size() / 2);
    for (size_t i = 0; i < mesh.Vertices.size(); i++)
    {
      mesh.Vertices[i].x = jsonPoints[2*i];
      mesh.Vertices[i].y = jsonPoints[2*i + 1];
    }

    // Extract cell data
    const auto jsonCells = json["Cells"];
    mesh.Cells.resize(jsonCells.size() / 3);
    for (size_t i = 0; i < mesh.Cells.size(); i++)
    {
      mesh.Cells[i].v0 = jsonCells[3*i];
      mesh.Cells[i].v1 = jsonCells[3*i + 1];
      mesh.Cells[i].v2 = jsonCells[3*i + 2];
    }
  }

  // Write 2D mesh to JSON file
  static void Write(const Mesh2D& mesh, std::string fileName)
  {
    std::cout << "JSON: Writing 2D mesh to file " << fileName << std::endl;

    // Generate JSON data
    auto jsonPoints = nlohmann::json::array();
    auto jsonCells = nlohmann::json::array();
    for (const auto p: mesh.Vertices)
    {
      jsonPoints.push_back(p.x);
      jsonPoints.push_back(p.y);
    }
    for (const auto T: mesh.Cells)
    {
      jsonCells.push_back(T.v0);
      jsonCells.push_back(T.v1);
      jsonCells.push_back(T.v2);
    }

    // Set JSON data
    nlohmann::json json;
    json["Type"] = "Mesh2D";
    json["Points"] = jsonPoints;
    json["Cells"] = jsonCells;

    // Write to file
    std::ofstream f(fileName);
    f << json;
  }

  // Read 3D mesh from JSON file
  static void Read(Mesh3D &mesh, std::string fileName)
  {
    std::cout << "JSON: Reading 3D mesh from file " << fileName << std::endl;

    // Read data from file
    std::ifstream f(fileName);
    nlohmann::json json;
    f >> json;

    // Check file type
    if (json.find("Type") == json.end() || json["Type"] != "Mesh3D")
      throw std::runtime_error("Not a DTCC Mesh3D JSON file.");

    // Extract point data
    const auto jsonPoints = json["Points"];
    mesh.Vertices.resize(jsonPoints.size() / 3);
    for (size_t i = 0; i < mesh.Vertices.size(); i++)
    {
      mesh.Vertices[i].x = jsonPoints[3*i];
      mesh.Vertices[i].y = jsonPoints[3*i + 1];
      mesh.Vertices[i].z = jsonPoints[3*i + 2];
    }

    // Extract cell data
    const auto jsonCells = json["Cells"];
    mesh.Cells.resize(jsonCells.size() / 4);
    for (size_t i = 0; i < mesh.Cells.size(); i++)
    {
      mesh.Cells[i].v0 = jsonCells[4*i];
      mesh.Cells[i].v1 = jsonCells[4*i + 1];
      mesh.Cells[i].v2 = jsonCells[4*i + 2];
      mesh.Cells[i].v3 = jsonCells[4*i + 3];
    }
  }

  // Write 3D mesh to JSON file
  static void Write(const Mesh3D& mesh, std::string fileName)
  {
    std::cout << "JSON: Writing 3D mesh to file " << fileName << std::endl;

    // Generate JSON data
    auto jsonPoints = nlohmann::json::array();
    auto jsonCells = nlohmann::json::array();
    for (const auto p: mesh.Vertices)
    {
      jsonPoints.push_back(p.x);
      jsonPoints.push_back(p.y);
      jsonPoints.push_back(p.z);
    }
    for (const auto T: mesh.Cells)
    {
      jsonCells.push_back(T.v0);
      jsonCells.push_back(T.v1);
      jsonCells.push_back(T.v2);
      jsonCells.push_back(T.v3);
    }

    // Set JSON data
    nlohmann::json json;
    json["Type"] = "Mesh3D";
    json["Points"] = jsonPoints;
    json["Cells"] = jsonCells;

    // Write to file
    std::ofstream f(fileName);
    f << json;
  }

  // Read 3D surface from JSON file
  static void Read(Surface3D &surface, std::string fileName)
  {
    std::cout << "JSON: Reading 3D surface from file " << fileName << std::endl;

    // Read data from file
    std::ifstream f(fileName);
    nlohmann::json json;
    f >> json;

    // Check file type
    if (json.find("Type") == json.end() || json["Type"] != "Surface3D")
      throw std::runtime_error("Not a DTCC Surface3D JSON file.");

    // Extract point data
    const auto jsonPoints = json["Points"];
    surface.Vertices.resize(jsonPoints.size() / 3);
    for (size_t i = 0; i < surface.Vertices.size(); i++)
    {
      surface.Vertices[i].x = jsonPoints[3*i];
      surface.Vertices[i].y = jsonPoints[3*i + 1];
      surface.Vertices[i].z = jsonPoints[3*i + 2];
    }

    // Extract cell data
    const auto jsonCells = json["Cells"];
    surface.Cells.resize(jsonCells.size() / 3);
    for (size_t i = 0; i < surface.Cells.size(); i++)
    {
      surface.Cells[i].v0 = jsonCells[3*i];
      surface.Cells[i].v1 = jsonCells[3*i + 1];
      surface.Cells[i].v2 = jsonCells[3*i + 2];
    }
  }

  // Write 3D surface to JSON file
  static void Write(const Surface3D& surface, std::string fileName)
  {
    std::cout << "JSON: Writing 3D surface to file " << fileName << std::endl;

    // Generate JSON data
    auto jsonPoints = nlohmann::json::array();
    auto jsonCells = nlohmann::json::array();
    for (const auto p: surface.Vertices)
    {
      jsonPoints.push_back(p.x);
      jsonPoints.push_back(p.y);
      jsonPoints.push_back(p.z);
    }
    for (const auto T: surface.Cells)
    {
      jsonCells.push_back(T.v0);
      jsonCells.push_back(T.v1);
      jsonCells.push_back(T.v2);
    }

    // Set JSON data
    nlohmann::json json;
    json["Type"] = "Surface3D";
    json["Points"] = jsonPoints;
    json["Cells"] = jsonCells;

    // Write to file
    std::ofstream f(fileName);
    f << json;
  }

private:
  // Try reading string value
  static std::string TryReadString(std::string key, const nlohmann::json &json)
  {
    if (json.find(key) == json.end())
      throw std::runtime_error("Missing field '" + key + "' in JSON file.");
    return json[key];
  }

  // Try reading double value
  static double TryReadDouble(std::string key, const nlohmann::json &json)
  {
    if (json.find(key) == json.end())
      throw std::runtime_error("Missing field '" + key + "' in JSON file.");
    return json[key];
  }

  // Try reading double value
  static bool TryReadBool(std::string key, const nlohmann::json &json)
  {
    if (json.find(key) == json.end())
      throw std::runtime_error("Missing field '" + key + "' in JSON file.");
    return json[key];
  }

  // Try reading integer value
  static int TryReadInt(std::string key, const nlohmann::json &json)
  {
    if (json.find(key) == json.end())
      throw std::runtime_error("Missing field '" + key + "' in JSON file.");
    return json[key];
  }
};

} // namespace DTCC

#endif
