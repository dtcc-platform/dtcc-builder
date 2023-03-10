// Copyright (C) 2020-2021 Anders Logg, Anton J Olsson
// Licensed under the MIT License

#ifndef DTCC_JSON_H
#define DTCC_JSON_H

#include <datamodel/CityModel.h>
#include <datamodel/District.h>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <utility>

#include "BoundingBox.h"
#include "Color.h"
#include "ColorMap.h"
#include "Grid.h"
#include "GridField.h"
#include "GridVectorField.h"
#include "Mesh.h"
#include "Parameters.h"
#include "PointCloud.h"
#include "Surface.h"
#include "cityjson/CityJSON.h"
#include "datamodel/RoadNetwork.h"

namespace DTCC_BUILDER
{

  class JSON
  {
  public:

    // Utility functions for JSON input/output
    #include "JSONUtils.h"

    /// Deserialize Parameters
    static void Deserialize(Parameters& parameters, const nlohmann::json& json)
    {
      CheckType("Parameters", json);

      // Read parameter values
      for (auto &it : parameters.Map)
      {
        if (!HasKey(it.first, json))
          warning("Missing parameter \"" + it.first +
                  "\", using default value");
        else
        {
          switch (it.second.Type)
          {
          case ParameterType::Bool:
            it.second = static_cast<bool>(json[it.first]);
            break;
          case ParameterType::Int:
            it.second = static_cast<int>(json[it.first]);
            break;
          case ParameterType::Float:
            it.second = static_cast<double>(json[it.first]);
            break;
          case ParameterType::String:
            std::string valueString = json[it.first];
            it.second = valueString;
            break;
          }
        }
      }

      // Report unknown parameters
      for (auto it = json.begin(); it != json.end(); ++it)
      {
        if (!parameters.HasKey(it.key()) && it.key() != "Type")
          warning("Unknown parameter \"" + it.key() + "\", ignoring");
      }
    };

    /// Serialize Parameters
    static void Serialize(const Parameters& parameters, nlohmann::json& json)
    {
      json["Type"] = "Parameters";
      for (auto &it : parameters.Map)
      {
        switch (it.second.Type)
        {
        case ParameterType::Bool:
          json[it.first] = static_cast<bool>(it.second);
          break;
        case ParameterType::Int:
          json[it.first] = static_cast<int>(it.second);
          break;
        case ParameterType::Float:
          json[it.first] = static_cast<double>(it.second);
          break;
        case ParameterType::String:
          json[it.first] = static_cast<std::string>(it.second);
          break;
        }
      }
    }

    /// Serialize timings (special-case, only output supported)
    static void
    Serialize(std::map<std::string, std::pair<double, size_t>> timings,
              nlohmann::json &json)
    {
      for (const auto &it : timings)
      {
        const std::string task = it.first;
        const double total = it.second.first;
        const size_t count = it.second.second;
        const double mean = total / static_cast<double>(count);
        auto jsonTiming = nlohmann::json::object();
        jsonTiming["Total"] = total;
        jsonTiming["Count"] = count;
        jsonTiming["Mean"] = mean;
        json[task] = jsonTiming;
      }
    }

    /// Deserialize BoundingBox2D
    static void Deserialize(BoundingBox2D& boundingBox, const nlohmann::json& json)
    {
      CheckType("BoundingBox2D", json);
      boundingBox.P.x = ToDouble("x", json["P"]);
      boundingBox.P.y = ToDouble("y", json["P"]);
      boundingBox.Q.x = ToDouble("x", json["Q"]);
      boundingBox.Q.y = ToDouble("y", json["Q"]);
    }

    /// Serialize BoundingBox2D
    static void Serialize(const BoundingBox2D& boundingBox, nlohmann::json& json)
    {
      auto jsonP = nlohmann::json::object();
      jsonP["x"] = boundingBox.P.x;
      jsonP["y"] = boundingBox.P.y;
      auto jsonQ = nlohmann::json::object();
      jsonQ["x"] = boundingBox.Q.x;
      jsonQ["y"] = boundingBox.Q.y;
      json["Type"] = "BoundingBox2D";
      json["P"] = jsonP;
      json["Q"] = jsonQ;
    }

    /// Deserialize BoundingBox3D
    static void Deserialize(BoundingBox3D& boundingBox, const nlohmann::json& json)
    {
      CheckType("BoundingBox3D", json);
      boundingBox.P.x = ToDouble("x", json["P"]);
      boundingBox.P.y = ToDouble("y", json["P"]);
      boundingBox.P.z = ToDouble("z", json["P"]);
      boundingBox.Q.x = ToDouble("x", json["Q"]);
      boundingBox.Q.y = ToDouble("y", json["Q"]);
      boundingBox.Q.z = ToDouble("z", json["Q"]);
    }

    /// Serialize BoundingBox3D
    static void Serialize(const BoundingBox3D& boundingBox, nlohmann::json& json)
    {
      auto jsonP = nlohmann::json::object();
      jsonP["x"] = boundingBox.P.x;
      jsonP["y"] = boundingBox.P.y;
      jsonP["z"] = boundingBox.P.z;
      auto jsonQ = nlohmann::json::object();
      jsonQ["x"] = boundingBox.Q.x;
      jsonQ["y"] = boundingBox.Q.y;
      jsonQ["z"] = boundingBox.Q.z;
      json["Type"] = "BoundingBox3D";
      json["P"] = jsonP;
      json["Q"] = jsonQ;
    }

    /// Deserialize PointCloud (not including colors and classification)
    static void Deserialize(PointCloud &pointCloud, const nlohmann::json &json)
    {
      CheckType("PointCloud", json);
      Deserialize(pointCloud.BoundingBox, json["BoundingBox"]);
      const auto jsonPoints = json["Points"];
      pointCloud.Points.resize(jsonPoints.size() / 3);
      for (size_t i = 0; i < pointCloud.Points.size(); i++)
      {
        pointCloud.Points[i].x = jsonPoints[3 * i];
        pointCloud.Points[i].y = jsonPoints[3 * i + 1];
        pointCloud.Points[i].z = jsonPoints[3 * i + 2];
      }
    }

    /// Serialize PointCloud
    static void Serialize(const PointCloud &pointCloud, nlohmann::json &json)
    {
      json["Type"] = "PointCloud";
      auto jsonPoints = nlohmann::json::array();
      for (const auto p : pointCloud.Points)
      {
        jsonPoints.push_back(p.x);
        jsonPoints.push_back(p.y);
        jsonPoints.push_back(p.z);
      }
      json["Points"] = jsonPoints;
    }

    /// Deserialize Grid2D
    static void Deserialize(Grid2D& grid, const nlohmann::json& json)
    {
      CheckType("Grid2D", json);
      Deserialize(grid.BoundingBox, json["BoundingBox"]);
      grid.XSize = ToUnsignedInt("XSize", json);
      grid.YSize = ToUnsignedInt("YSize", json);
      grid.XStep = ToDouble("XStep", json);
      grid.YStep = ToDouble("YStep", json);
    }

    /// Serialize Grid2D
    static void Serialize(const Grid2D& grid, nlohmann::json& json)
    {
      auto jsonBoundingBox = nlohmann::json::object();
      Serialize(grid.BoundingBox, jsonBoundingBox);
      json["Type"] = "Grid2D";
      json["BoundingBox"] = jsonBoundingBox;
      json["XSize"] = grid.XSize;
      json["YSize"] = grid.YSize;
      json["XStep"] = grid.XStep;
      json["YStep"] = grid.YStep;
    }

    /// Deserialize Grid3D
    static void Deserialize(Grid3D& grid, const nlohmann::json& json)
    {
      CheckType("Grid3D", json);
      Deserialize(grid.BoundingBox, json["BoundingBox"]);
      grid.XSize = ToUnsignedInt("XSize", json);
      grid.YSize = ToUnsignedInt("YSize", json);
      grid.ZSize = ToUnsignedInt("ZSize", json);
      grid.XStep = ToDouble("XStep", json);
      grid.YStep = ToDouble("YStep", json);
      grid.ZStep = ToDouble("ZStep", json);
    }

    /// Serialize Grid3D
    static void Serialize(const Grid3D& grid, nlohmann::json& json)
    {
      auto jsonBoundingBox = nlohmann::json::object();
      Serialize(grid.BoundingBox, jsonBoundingBox);
      json["Type"] = "Grid3D";
      json["BoundingBox"] = jsonBoundingBox;
      json["XSize"] = grid.XSize;
      json["YSize"] = grid.YSize;
      json["ZSize"] = grid.ZSize;
      json["XStep"] = grid.XStep;
      json["YStep"] = grid.YStep;
      json["ZStep"] = grid.ZStep;
    }

    /// Deserialize Mesh2D
    static void Deserialize(Mesh2D& mesh, const nlohmann::json& json)
    {
      CheckType("Mesh2D", json);
      const auto jsonVertices = json["Vertices"];
      mesh.Vertices.resize(jsonVertices.size() / 2);
      for (size_t i = 0; i < mesh.Vertices.size(); i++)
      {
        mesh.Vertices[i].x = jsonVertices[2*i];
        mesh.Vertices[i].y = jsonVertices[2*i + 1];
      }
      const auto jsonCells = json["Cells"];
      mesh.Cells.resize(jsonCells.size() / 3);
      for (size_t i = 0; i < mesh.Cells.size(); i++)
      {
        mesh.Cells[i].v0 = jsonCells[3*i];
        mesh.Cells[i].v1 = jsonCells[3*i + 1];
        mesh.Cells[i].v2 = jsonCells[3*i + 2];
      }
      const auto jsonMarkers = json["Markers"];
      mesh.Markers.resize(jsonMarkers.size());
      for (size_t i = 0; i < mesh.Markers.size(); i++)
        mesh.Markers[i] = jsonMarkers[i];
    }

    /// Serialize Mesh2D and its offset/origin
    static void
    Serialize(const Mesh2D &mesh, nlohmann::json &json, const Point2D &origin)
    {
      auto jsonVertices = nlohmann::json::array();
      for (const auto p: mesh.Vertices)
      {
        jsonVertices.push_back(p.x);
        jsonVertices.push_back(p.y);
      }
      auto jsonCells = nlohmann::json::array();
      for (const auto T: mesh.Cells)
      {
        jsonCells.push_back(T.v0);
        jsonCells.push_back(T.v1);
        jsonCells.push_back(T.v2);
      }
      json["Type"] = "Mesh2D";
      json["Vertices"] = jsonVertices;
      json["Cells"] = jsonCells;
      json["Markers"] = mesh.Markers;
      json["Origin"] = {{"x", origin.x}, {"y", origin.y}};
    }

    /// Deserialize Mesh3D
    static void Deserialize(Mesh3D& mesh, const nlohmann::json& json)
    {
      CheckType("Mesh3D", json);
      const auto jsonVertices = json["Vertices"];
      mesh.Vertices.resize(jsonVertices.size() / 3);
      for (size_t i = 0; i < mesh.Vertices.size(); i++)
      {
        mesh.Vertices[i].x = jsonVertices[3*i];
        mesh.Vertices[i].y = jsonVertices[3*i + 1];
        mesh.Vertices[i].z = jsonVertices[3*i + 2];
      }
      const auto jsonCells = json["Cells"];
      mesh.Cells.resize(jsonCells.size() / 4);
      for (size_t i = 0; i < mesh.Cells.size(); i++)
      {
        mesh.Cells[i].v0 = jsonCells[4*i];
        mesh.Cells[i].v1 = jsonCells[4*i + 1];
        mesh.Cells[i].v2 = jsonCells[4*i + 2];
        mesh.Cells[i].v3 = jsonCells[4*i + 3];
      }
      const auto jsonMarkers = json["Markers"];
      mesh.Markers.resize(jsonMarkers.size());
      for (size_t i = 0; i < mesh.Markers.size(); i++)
        mesh.Markers[i] = jsonMarkers[i];
    }

    /// Serialize Mesh3D and its offset/origin
    static void
    Serialize(const Mesh3D &mesh, nlohmann::json &json, const Point2D &origin)
    {
      auto jsonVertices = nlohmann::json::array();
      for (const auto p: mesh.Vertices)
      {
        jsonVertices.push_back(p.x);
        jsonVertices.push_back(p.y);
        jsonVertices.push_back(p.z);
      }
      auto jsonCells = nlohmann::json::array();
      for (const auto T: mesh.Cells)
      {
        jsonCells.push_back(T.v0);
        jsonCells.push_back(T.v1);
        jsonCells.push_back(T.v2);
        jsonCells.push_back(T.v3);
      }
      json["Type"] = "Mesh3D";
      json["Vertices"] = jsonVertices;
      json["Cells"] = jsonCells;
      json["Markers"] = mesh.Markers;
      json["Origin"] = {{"x", origin.x}, {"y", origin.y}};
    }

    /// Deserialize Surface2D
    static void Deserialize(Surface2D &surface, const nlohmann::json& json)
    {
      CheckType("Surface3D", json);
      const auto jsonVertices = json["Vertices"];
      surface.Vertices.resize(jsonVertices.size() / 2);
      for (size_t i = 0; i < surface.Vertices.size(); i++)
      {
        surface.Vertices[i].x = jsonVertices[2*i];
        surface.Vertices[i].y = jsonVertices[2*i + 1];
      }
      const auto jsonEdges = json["Edges"];
      surface.Edges.resize(jsonEdges.size() / 2);
      for (size_t i = 0; i < surface.Edges.size(); i++)
      {
        surface.Edges[i].v0 = jsonEdges[2 * i];
        surface.Edges[i].v1 = jsonEdges[2 * i + 1];
      }
      const auto jsonNormals = json["Normals"];
      surface.Normals.resize(jsonNormals.size() / 2);
      for (size_t i = 0; i < surface.Normals.size(); i++)
      {
        surface.Normals[i].x = jsonNormals[2 * i];
        surface.Normals[i].y = jsonNormals[2 * i + 1];
      }
    }

    /// Serialize Surface2D and its offset/origin
    static void Serialize(const Surface2D &surface,
                          nlohmann::json &json,
                          const Point2D &origin)
    {
      auto jsonVertices = nlohmann::json::array();
      for (const auto p: surface.Vertices)
      {
        jsonVertices.push_back(p.x);
        jsonVertices.push_back(p.y);
      }
      auto jsonEdges = nlohmann::json::array();
      for (const auto T : surface.Edges)
      {
        jsonEdges.push_back(T.v0);
        jsonEdges.push_back(T.v1);
      }
      auto jsonNormals = nlohmann::json::array();
      for (const auto p : surface.Normals)
      {
        jsonNormals.push_back(p.x);
        jsonNormals.push_back(p.y);
      }
      json["Type"] = "Surface2D";
      json["Vertices"] = jsonVertices;
      json["Edges"] = jsonEdges;
      json["Normals"] = jsonNormals;
      json["Origin"] = {{"x", origin.x}, {"y", origin.y}};
    }

    /// Deserialize Surface3D
    static void Deserialize(Surface3D &surface, const nlohmann::json& json)
    {
      CheckType("Surface3D", json);
      const auto jsonVertices = json["Vertices"];
      surface.Vertices.resize(jsonVertices.size() / 3);
      for (size_t i = 0; i < surface.Vertices.size(); i++)
      {
        surface.Vertices[i].x = jsonVertices[3*i];
        surface.Vertices[i].y = jsonVertices[3*i + 1];
        surface.Vertices[i].z = jsonVertices[3*i + 2];
      }
      const auto jsonFaces = json["Faces"];
      surface.Faces.resize(jsonFaces.size() / 3);
      for (size_t i = 0; i < surface.Faces.size(); i++)
      {
        surface.Faces[i].v0 = jsonFaces[3 * i];
        surface.Faces[i].v1 = jsonFaces[3 * i + 1];
        surface.Faces[i].v2 = jsonFaces[3 * i + 2];
      }
      const auto jsonNormals = json["Normals"];
      surface.Normals.resize(jsonNormals.size() / 3);
      for (size_t i = 0; i < surface.Normals.size(); i++)
      {
        surface.Normals[i].x = jsonNormals[3 * i];
        surface.Normals[i].y = jsonNormals[3 * i + 1];
        surface.Normals[i].z = jsonNormals[3 * i + 2];
      }
    }

    /// Serialize Surface3D and its offset/origin
    static void Serialize(const Surface3D &surface,
                          nlohmann::json &json,
                          const Point2D &origin)
    {
      auto jsonVertices = nlohmann::json::array();
      for (const auto p: surface.Vertices)
      {
        jsonVertices.push_back(p.x);
        jsonVertices.push_back(p.y);
        jsonVertices.push_back(p.z);
      }
      auto jsonFaces = nlohmann::json::array();
      for (const auto T : surface.Faces)
      {
        jsonFaces.push_back(T.v0);
        jsonFaces.push_back(T.v1);
        jsonFaces.push_back(T.v2);
      }
      auto jsonNormals = nlohmann::json::array();
      for (const auto p : surface.Normals)
      {
        jsonNormals.push_back(p.x);
        jsonNormals.push_back(p.y);
        jsonNormals.push_back(p.z);
      }
      json["Type"] = "Surface3D";
      json["Vertices"] = jsonVertices;
      json["Faces"] = jsonFaces;
      json["Normals"] = jsonNormals;
      json["Origin"] = {{"x", origin.x}, {"y", origin.y}};
    }

    /// Deserialize GridField2D
    static void Deserialize(GridField2D& field, const nlohmann::json& json)
    {
      CheckType("GridField2D", json);
      Deserialize(field.Grid, json["Grid"]);
      const auto jsonValues = json["Values"];
      field.Values.resize(jsonValues.size());
      for (size_t i = 0; i < field.Values.size(); i++)
        field.Values[i] = jsonValues[i];
    }

    /// Serialize GridField2D and its offset/origin
    static void Serialize(const GridField2D &field,
                          nlohmann::json &json,
                          const Point2D &origin)
    {
      auto jsonGrid = nlohmann::json::object();
      Serialize(field.Grid, jsonGrid);
      json["Type"] = "GridField2D";
      json["Grid"] = jsonGrid;
      json["Values"] = field.Values;
      json["Origin"] = {{"x", origin.x}, {"y", origin.y}};
    }

    /// Deserialize GridField3D
    static void Deserialize(GridField3D& field, const nlohmann::json& json)
    {
      CheckType("GridField3D", json);
      Deserialize(field.Grid, json["Grid"]);
      const auto jsonValues = json["Values"];
      field.Values.resize(jsonValues.size());
      for (size_t i = 0; i < field.Values.size(); i++)
        field.Values[i] = jsonValues[i];
    }

    /// Serialize GridField3D and its offset/origin
    static void Serialize(const GridField3D &field,
                          nlohmann::json &json,
                          const Point2D &origin)
    {
      auto jsonGrid = nlohmann::json::object();
      Serialize(field.Grid, jsonGrid);
      json["Type"] = "GridField3D";
      json["Grid"] = jsonGrid;
      json["Values"] = field.Values;
      json["Origin"] = {{"x", origin.x}, {"y", origin.y}};
    }

    /// Deserialize GridVectorField2D
    static void Deserialize(GridVectorField2D& field, const nlohmann::json& json)
    {
      CheckType("GridVectorField2D", json);
      Deserialize(field.Grid, json["Grid"]);
      const auto jsonValues = json["Values"];
      field.Values.resize(jsonValues.size());
      for (size_t i = 0; i < field.Values.size(); i++)
        field.Values[i] = jsonValues[i];
    }

    /// Serialize GridVectorField2D and its offset/origin
    static void Serialize(const GridVectorField2D &field,
                          nlohmann::json &json,
                          const Point2D &origin)
    {
      auto jsonGrid = nlohmann::json::object();
      Serialize(field.Grid, jsonGrid);
      json["Type"] = "GridVectorField2D";
      json["Grid"] = jsonGrid;
      json["Values"] = field.Values;
      json["Origin"] = {{"x", origin.x}, {"y", origin.y}};
    }

    /// Deserialize GridVectorField3D
    static void Deserialize(GridVectorField3D& field, const nlohmann::json& json)
    {
      CheckType("GridVectorField3D", json);
      Deserialize(field.Grid, json["Grid"]);
      const auto jsonValues = json["Values"];
      field.Values.resize(jsonValues.size());
      for (size_t i = 0; i < field.Values.size(); i++)
        field.Values[i] = jsonValues[i];
    }

    /// Serialize GridVectorField3D and its offset/origin
    static void Serialize(const GridVectorField3D &field,
                          nlohmann::json &json,
                          const Point2D &origin)
    {
      auto jsonGrid = nlohmann::json::object();
      Serialize(field.Grid, jsonGrid);
      json["Type"] = "GridVectorField3D";
      json["Grid"] = jsonGrid;
      json["Values"] = field.Values;
      json["Origin"] = {{"x", origin.x}, {"y", origin.y}};
    }

    // FIXME: Add separate serialize/deserialize for Building and use here.

    /// Deserialize CityModel
    static void Deserialize(CityModel &cityModel, const nlohmann::json& json)
    {
      CheckType("CityModel", json);
      cityModel.Name = json["Name"];
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
        auto jsonRoofPoints = jsonBuilding["RoofPoints"];
        cityModel.Buildings[i].RoofPoints.resize(jsonRoofPoints.size());
        for (size_t j = 0; j < jsonRoofPoints.size(); j++)
        {
          cityModel.Buildings[i].RoofPoints[j].x = jsonRoofPoints[j]["x"];
          cityModel.Buildings[i].RoofPoints[j].y = jsonRoofPoints[j]["y"];
          cityModel.Buildings[i].RoofPoints[j].z = jsonRoofPoints[j]["z"];
        }
        auto jsonGroundPoints = jsonBuilding["GroundPoints"];
        cityModel.Buildings[i].GroundPoints.resize(jsonGroundPoints.size());
        for (size_t j = 0; j < jsonGroundPoints.size(); j++)
        {
          cityModel.Buildings[i].GroundPoints[j].x = jsonGroundPoints[j]["x"];
          cityModel.Buildings[i].GroundPoints[j].y = jsonGroundPoints[j]["y"];
          cityModel.Buildings[i].GroundPoints[j].z = jsonGroundPoints[j]["z"];
        }
        cityModel.Buildings[i].Height = jsonBuilding["Height"];
        cityModel.Buildings[i].GroundHeight = jsonBuilding["GroundHeight"];
        cityModel.Buildings[i].UUID = jsonBuilding["UUID"];
        cityModel.Buildings[i].SHPFileID = jsonBuilding["SHPFileID"];
        cityModel.Buildings[i].error = jsonBuilding["Error"];
      }
    }

    /// Serialize CityModel and its offset/origin
    static void Serialize(const CityModel &cityModel,
                          nlohmann::json &json,
                          const Point2D &origin)
    {
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
        jsonBuilding["RoofPoints"] = nlohmann::json::array();
        for (auto const &point : building.RoofPoints)
        {
          auto jsonPoint = nlohmann::json::object();
          jsonPoint["x"] = point.x;
          jsonPoint["y"] = point.y;
          jsonPoint["z"] = point.z;
          jsonBuilding["RoofPoints"].push_back(jsonPoint);
        }
        jsonBuilding["GroundPoints"] = nlohmann::json::array();
        for (auto const &point : building.GroundPoints)
        {
          auto jsonPoint = nlohmann::json::object();
          jsonPoint["x"] = point.x;
          jsonPoint["y"] = point.y;
          jsonPoint["z"] = point.z;
          jsonBuilding["GroundPoints"].push_back(jsonPoint);
        }
        jsonBuilding["Height"] = building.Height;
        jsonBuilding["GroundHeight"] = building.GroundHeight;
        jsonBuilding["UUID"] = building.UUID;
        // Uncomment for debugging
        // jsonBuilding["debugID"] = building.debugID;
        jsonBuilding["SHPFileID"] = building.SHPFileID;
        jsonBuilding["Error"] = building.error;
        jsonBuildings.push_back(jsonBuilding);
      }
      json["Type"] = "CityModel";
      json["Name"] = cityModel.Name;
      json["Buildings"] = jsonBuildings;
      json["Origin"] = {{"x", origin.x}, {"y", origin.y}};
    }

    /// Serialize CityJSON
    static void Serialize(const CityJSON &cityJson, nlohmann::json &json)
    {
      // Serializing city objects
      auto jsonBuilding = nlohmann::json::object();
      for(auto const& cityObject : cityJson.CityObjects)
      {
        //jsonBuilding[cityObject.ID]
        // Storing the id of the object
        std::string id = cityObject.ID;
        auto objID = nlohmann::json::object();
        json["CityObjects"][id]=objID;

        // Storing the attributes object
        auto attributesObj = nlohmann::json::object();
        attributesObj["measuredHeight"]=cityObject.ObjectAttributes.MeasuredHeight;
        // TODO: Extend this with additional attributes
        json["CityObjects"][id]["attributes"]=attributesObj;

        // Storing the geometry object
        auto geometryObj = nlohmann::json::object();
        geometryObj["lod"]=cityObject.ObjectGeometry.LOD;
        geometryObj["type"]=CityObject::Geometry::GeometryTypeToString(cityObject.ObjectGeometry.Type);
        // TODO: Extend this with additional geometry fields

        //json["CityObjects"][id]["geometry"]=geometryObj;

        auto geometryArray=nlohmann::json::array();

        if(cityObject.ObjectGeometry.Type==CityObject::Geometry::GeometryType::Solid)
        {
          // Storing the boundaries array for solid meshes
          auto boundariesExternalArray = nlohmann::json::array();
          auto boundariesInternalArray = nlohmann::json::array();

          for (auto const &boundary : cityObject.ObjectGeometry.Boundaries)
          {
            auto singleBoundaryArray = nlohmann::json::array();
            auto boundariesIDs = nlohmann::json::array();
            for (auto const &id : boundary.BoundariesIDs)
            {
              boundariesIDs.push_back(id);
            }
            singleBoundaryArray.push_back(boundariesIDs);
            boundariesInternalArray.push_back(singleBoundaryArray);
          }

          boundariesExternalArray.push_back(boundariesInternalArray);

          geometryObj["boundaries"] = boundariesExternalArray;
          geometryArray.push_back(geometryObj);
          json["CityObjects"][id]["geometry"] = geometryArray;

          // TODO: Change this accordingly to support different types
          objID["type"] = "Building";
        }
        else //Treat everything as multisurface for now
        {
            auto externalArray = nlohmann::json::array();

            for(auto const& boundary : cityObject.ObjectGeometry.Boundaries)
            {
              auto boundariesInternalArray = nlohmann::json::array();
              auto boundariesArray = nlohmann::json::array();
              for(auto const& id: boundary.BoundariesIDs)
              {
                boundariesArray.push_back(id);
              }
              boundariesInternalArray.push_back(boundariesArray);
              externalArray.push_back(boundariesInternalArray);
            }

            geometryObj["boundaries"] = externalArray;
            geometryArray.push_back(geometryObj);
            json["CityObjects"][id]["geometry"] = geometryArray;
        }
      }

      //Serializing vertices
      auto jsonVertices= nlohmann::json::array();
      for(auto const& vertex: cityJson.Vertices)
      {
        auto jsonVertex = nlohmann::json::array();

        jsonVertex.push_back(vertex.x);
        jsonVertex.push_back(vertex.y);
        jsonVertex.push_back(vertex.z);
        jsonVertices.push_back(jsonVertex);
      }
      json["vertices"] = jsonVertices;




      json["type"] = "CityJSON";
      json["version"] = "1.0";
    }

    /// Deserialize CityJSON
    static void Deserialize(CityJSON &cityJson, const nlohmann::json &json)
    {
      //CheckType("CityJSON",json);
      //Getting vertices
      uint verticesNum = json["vertices"].size();
      for (uint i = 0; i < verticesNum; i++)
      {
        // debug
        //std::cout << json["vertices"][i] << std::endl;
        auto vertex = json["vertices"][i];

        assert(vertex.size() == 3 &&
               "Attempted to read non 3D vertex from json file");

        cityJson.Vertices.push_back(Point3D(json["vertices"][i][0],
                                   json["vertices"][i][1],
                                   json["vertices"][i][2]));
      }

      //Getting city objects
      auto cityObjectsField = json["CityObjects"].get<nlohmann::json::object_t>();
      for (auto cityObjectIterator = cityObjectsField.begin();
           cityObjectIterator != cityObjectsField.end(); ++cityObjectIterator)
      {
        //debug
        //std::cout << cityObjectIterator->first << std::endl; // getting City Object ID here

        // Creating a new object in stack, we're going to store the information
        // here and later store this into cityObjects vector
        CityObject NewObj = CityObject(cityObjectIterator->first);

        auto jsonCityObj = json["CityObjects"][cityObjectIterator->first];
        //std::cout<<jsonCityObj<<std::endl;
        // Getting attributes -- TODO Make this error free
        CityObject::Attributes attributes = CityObject::Attributes();
        auto attributesField = json.find("attributes");
        if(attributesField!=json.end())
        {
          // TODO: Get more fields for attributes here
          //attributes = attributesField["measuredHeight"];
          attributes.MeasuredHeight = jsonCityObj["attributes"]["measuredHeight"];
        }

        NewObj.ObjectAttributes=attributes;

        // Storing geometry
        CityObject::Geometry geometry = CityObject::Geometry();
        auto jsonGeometryArray = json["CityObjects"][NewObj.ID]["geometry"][0];
        //std::cout << jsonGeometryArray << std::endl;

        // Storing lod settings
        NewObj.ObjectGeometry.LOD = jsonGeometryArray["lod"];
        std::cout << "LOD:"<<geometry.LOD<<std::endl;

        // Storing building type
        std::string tempGeomType = jsonGeometryArray["type"];
        std::cout<<tempGeomType<<std::endl;
        // TODO: correlate this to string-enum
        if(tempGeomType=="Solid")
        {
          NewObj.ObjectGeometry.Type = CityObject::Geometry::Solid;
        }
        else if (tempGeomType == "MultiSurface")
        {
          NewObj.ObjectGeometry.Type = CityObject::Geometry::MultiSurface;
          // TODO: this should be replaced by built-in functionality in struct
        }

        // Storing boundaries
        uint boundariesArraySize = jsonGeometryArray["boundaries"][0].size();

        if (NewObj.ObjectGeometry.Type == CityObject::Geometry::Solid)
        {
          for (uint i = 0; i < boundariesArraySize; i++)
          {
            auto boundariesJson = jsonGeometryArray["boundaries"][0][i];
            for (auto boundaryIter = boundariesJson.begin();
                 boundaryIter != boundariesJson.end(); ++boundaryIter)
            {
              // std::cout << "boundary iter:" << *boundaryIter << std::endl;
              auto internalJsonArray = boundariesJson[0];

              CityObject::Geometry::Boundary newBoundary;

              for (auto it = internalJsonArray.begin();
                   it != internalJsonArray.end(); ++it)
              {
                // std::cout << *it << std::endl;
                newBoundary.BoundariesIDs.push_back(*it);
              }

              NewObj.ObjectGeometry.Boundaries.push_back(newBoundary);
            }
          }
        }
        else //Assume everything else is multi-surface for now
        {
          for(auto& outer: jsonGeometryArray["boundaries"])
          {
            CityObject::Geometry::Boundary newBoundary;
            for(auto& inner:outer)
            {
                for (auto& id : inner)
                {
                    newBoundary.BoundariesIDs.push_back(id);
                }
            }
            NewObj.ObjectGeometry.Boundaries.push_back(newBoundary);
          }
        }

        // Getting type
        // TODO: move this whole impl to cityobj class function
        auto objectType = jsonCityObj.find("type");
        if(objectType!=jsonCityObj.end())
        {
          NewObj.ObjectType = CityObject::Building;
        }
        cityJson.CityObjects.push_back(NewObj);
      }

    }

    static void Serialize(const ColorMap &colorMap, nlohmann::json &json)
    {
      auto jsonColorMap = nlohmann::json::array();
      json["Type"] = "ColorMap";
      if (colorMap.mapType == Linear)
        json["colorMapType"] = "Linear";
      else
        json["colorMapType"] = "Discrete";
      for (const auto c : colorMap.Colors)
      {
        auto jsonColorMapEntry = nlohmann::json::array();
        jsonColorMapEntry.push_back(c.first);
        Color cl = c.second;
        jsonColorMapEntry.push_back(cl.R);
        jsonColorMapEntry.push_back(cl.G);
        jsonColorMapEntry.push_back(cl.B);
        jsonColorMap.push_back(jsonColorMapEntry);
      }
      json["map"] = jsonColorMap;
    }

    static void Deserialize(ColorMap& colorMap, const nlohmann::json& json)
    {
      CheckType("ColorMap", json);
      if (json["colorMapType"] == "Linear")
        colorMap.mapType = Linear;
      else
        colorMap.mapType = Discrete;
      const auto colorMapEntries = json["map"];
      for (auto c : colorMapEntries)
      {
        colorMap.InsertColor((double)c[0],Color((double)c[1],(double)c[2],(double)c[3]));
      }

    }

    /// Deserialize RoadNetwork
    static void Deserialize(RoadNetwork &road, const nlohmann::json &json)
    {
      CheckType("RoadNetwork", json);
      auto jsonRoadNetwork = json["RoadNetwork"];

      // Read vertex positions
      const auto jsonVertices = jsonRoadNetwork["Vertices"];
      road.Vertices.resize(jsonVertices.size() / 2);
      for (size_t i = 0; i < road.Vertices.size(); i++)
      {
        road.Vertices[i].x = jsonVertices[i * 2];
        road.Vertices[i].y = jsonVertices[i * 2 + 1];
      }

      // Read edge indices
      const auto jsonEdges = jsonRoadNetwork["Edges"];
      road.Edges.resize(jsonEdges.size());
      for (size_t i = 0; i < road.Edges.size(); i++)
        road.Edges[i] = std::make_pair(jsonEdges[i][0], jsonEdges[i][1]);

      // Read additional vertex values
      const auto jsonVertexValues = jsonRoadNetwork["VertexValues"];
      for (auto it = jsonVertexValues.begin(); it != jsonVertexValues.end(); it++)
      {
        auto jsonVertexValuesArray = it.value();
        road.VertexValues.emplace(it.key(), std::vector<double>(jsonVertexValuesArray.size()));
        for (size_t i = 0; i < jsonVertexValuesArray.size(); i++)
        {
          road.VertexValues[it.key()][i] = jsonVertexValuesArray[i];
        }
      }

      // Read additional edge values
      const auto jsonEdgeValues = jsonRoadNetwork["EdgeValues"];
      for (auto it = jsonEdgeValues.begin(); it != jsonEdgeValues.end(); it++)
      {
        auto jsonEdgeValuesArray = it.value();
        road.EdgeValues.emplace(it.key(), std::vector<double>(jsonEdgeValuesArray.size()));
        for (size_t i = 0; i < jsonEdgeValuesArray.size(); i++)
        {
          road.EdgeValues[it.key()][i] = jsonEdgeValuesArray[i];
        }
      }

      // Read properties
      setJsonRoadProps(road.VertexProperties, jsonRoadNetwork,
                       "VertexProperties");
      setJsonRoadProps(road.EdgeProperties, jsonRoadNetwork, "EdgeProperties");
    }

    static void setJsonRoadProps(
        std::unordered_map<std::string, std::vector<std::string>> &propsMap,
        nlohmann::json &jsonRoadNetwork,
        const char *propsName)
    {
      const auto jsonProperties = jsonRoadNetwork[propsName];
      for (auto it = jsonProperties.begin(); it != jsonProperties.end(); it++)
      {
        auto jsonPropsArray = it.value();
        propsMap.emplace(it.key(),
                         jsonPropsArray.get<std::vector<std::string>>());
      }
    }

    /// Serialize RoadNetwork and its offset/origin
    static void
    Serialize(const RoadNetwork &road, nlohmann::json &json, Point2D origin)
    {
      auto jsonRoadNetwork = nlohmann::json::object();

      // Serialize Vertices
      auto jsonVertices = nlohmann::json::array();
      for (const auto p: road.Vertices)
      {
        jsonVertices.push_back(p.x);
        jsonVertices.push_back(p.y);
      }
      jsonRoadNetwork["Vertices"] = jsonVertices;

      // Serialize Edges
      auto jsonEdges = nlohmann::json::array();
      for (const auto &e : road.Edges)
      {
        nlohmann::json jsonEdge;
        jsonEdge.push_back(e.first);
        jsonEdge.push_back(e.second);
        jsonEdges.push_back(jsonEdge);
      }
      jsonRoadNetwork["Edges"] = jsonEdges;

      // Serialize VertexValues
      auto jsonVertexValues = nlohmann::json::object();
      for (auto it = road.VertexValues.begin(); it != road.VertexValues.end(); it++)
      {
        auto jsonVertexValuesArray = nlohmann::json::array();
        for (const auto v: it->second)
        {
          jsonVertexValuesArray.push_back(v);
        }
        jsonVertexValues[it->first] = jsonVertexValuesArray;
      }
      jsonRoadNetwork["VertexValues"] = jsonVertexValues;

      // Serialize EdgeValues
      auto jsonEdgeValues = nlohmann::json::object();
      for (auto it = road.EdgeValues.begin(); it != road.EdgeValues.end(); it++)
      {
        auto jsonEdgeValuesArray = nlohmann::json::array();
        for (const auto v: it->second)
        {
          jsonEdgeValuesArray.push_back(v);
        }
        jsonEdgeValues[it->first] = jsonEdgeValuesArray;
      }
      jsonRoadNetwork["EdgeValues"] = jsonEdgeValues;

      jsonRoadNetwork["VertexProperties"] =
          getJsonRoadProps(road.VertexProperties);
      jsonRoadNetwork["EdgeProperties"] = getJsonRoadProps(road.EdgeProperties);

      // Add to RoadNetwork object
      json["RoadNetwork"] = jsonRoadNetwork;
      // Add Type parameter
      json["Type"] = "RoadNetwork";
      // Add Origin parameter
      json["Origin"] = {{"x", origin.x}, {"y", origin.y}};
    }

    /// Deserialize Property.
    static void Deserialize(Property &property,
                            const nlohmann::json &json,
                            std::string uuid)
    {
      /// Get corresponding JSON object
      bool isProperty = json["Type"] == "Property";
      nlohmann::json jsonProperty =
          isProperty ? json
                     : GetObjectByAttribute("UUID", std::move(uuid),
                                            json["Properties"]);

      CheckType("Property", jsonProperty);

      DeserializeFootprint(property.Footprint, jsonProperty, "Footprint");
      property.FNR = jsonProperty["FNR"];
      property.UUID = jsonProperty["UUID"];
      /// Deserialize buildings.
      if (!property.Buildings.empty())
        return;
      auto jsonBuildings = jsonProperty["contains"];
      /// Deserialize properties.
      property.Buildings.resize(jsonBuildings.size());
      for (size_t i = 0; i < jsonBuildings.size(); ++i)
        Deserialize(property.Buildings[i], json, jsonBuildings[i]["UUID"]);
    }

    /// Serialize Property.
    static void Serialize(const Property &property,
                          nlohmann::json &json,
                          bool isTopObject = true,
                          bool serializeBuildings = true)
    {
      nlohmann::json jsonProperty;
      jsonProperty["UUID"] = property.UUID;
      jsonProperty["FNR"] = property.FNR;
      /// Serialize buildings.
      for (const auto &building : property.Buildings)
      {
        jsonProperty["contains"].push_back(building.UUID);
        if (serializeBuildings)
          Serialize(building, json, false);
      }
      SerializeArea(property, json, isTopObject, jsonProperty, "Property");
    }

    /// Serialize Building
    static void Serialize(const Building &building,
                          nlohmann::json &json,
                          bool isTopObject = true)
    {
      nlohmann::json jsonBuilding;
      jsonBuilding["UUID"] = building.UUID;
      jsonBuilding["BaseAreaID"] = building.BaseAreaID;
      jsonBuilding["PropertyFNR"] = building.PropertyFNR;
      jsonBuilding["Height"] = building.Height;
      jsonBuilding["GroundHeight"] = building.GroundHeight;
      jsonBuilding["PropertyUUID"] = building.PropertyUUID;
      jsonBuilding["Error"] = building.error;

      SerializeArea(building, json, isTopObject, jsonBuilding, "Building");
    }

    /// Deserialize Building.
    static void Deserialize(Building &building,
                            const nlohmann::json &json,
                            const std::string &uuid)
    {
      /// Get corresponding JSON object
      bool isBuilding = json["Type"] == "Building";
      nlohmann::json jsonBuilding =
          isBuilding ? json
                     : GetObjectByAttribute("UUID", uuid, json["Buildings"]);
      CheckType("Building", jsonBuilding);

      building.UUID = jsonBuilding["UUID"];
      building.PropertyUUID = jsonBuilding["PropertyUUID"];
      building.PropertyFNR = jsonBuilding["PropertyFNR"].get<int>();
      building.BaseAreaID = jsonBuilding["BaseAreaID"].get<std::string>();
      building.Height = jsonBuilding["Height"].get<double>();
      building.GroundHeight = jsonBuilding["GroundHeight"].get<double>();
      building.error = jsonBuilding["Error"];

      DeserializeFootprint(building.Footprint, jsonBuilding, "Footprint");
    }

    /// Serialize BaseArea.
    static void Serialize(const BaseArea &baseArea,
                          nlohmann::json &json,
                          bool isTopObject = true)
    {
      nlohmann::json jsonBaseArea;
      jsonBaseArea["area_id"] = baseArea.AreaID;
      jsonBaseArea["parent_id"] = baseArea.PrimaryAreaID;
      for (const auto &property : baseArea.Properties)
      {
        nlohmann::json jsonProperty;
        jsonProperty["UUID"] = property.UUID;
        jsonProperty["FNR"] = property.FNR;
        jsonBaseArea["contains"]["Properties"].push_back(jsonProperty);
        Serialize(property, json, false, false);
      }
      for (const auto &building : baseArea.Buildings)
      {
        jsonBaseArea["contains"]["Buildings"].push_back(building.UUID);
        Serialize(building, json, false);
      }
      SerializeArea(baseArea, json, isTopObject, jsonBaseArea, "BaseArea");
    }

    /// Deserialize BaseArea.
    static void Deserialize(BaseArea &baseArea,
                            const nlohmann::json &json,
                            std::string areaID)
    {
      /// Get corresponding JSON object
      bool isBaseArea = json["Type"] == "BaseArea";
      nlohmann::json jsonBaseArea =
          isBaseArea ? json
                     : GetObjectByAttribute("area_id", std::move(areaID),
                                            json["BaseAreas"]);
      CheckType("BaseArea", jsonBaseArea);

      baseArea.PrimaryAreaID = jsonBaseArea["parent_id"];
      DeserializeArea(baseArea, jsonBaseArea);
      nlohmann::json jsonBuildings = jsonBaseArea["contains"]["Buildings"];
      baseArea.Buildings.resize(jsonBuildings.size());
      for (size_t i = 0; i < jsonBuildings.size(); ++i)
        Deserialize(baseArea.Buildings[i], json,
                    jsonBuildings[i].get<std::string>());
      nlohmann::json jsonProperties = jsonBaseArea["contains"]["Properties"];
      baseArea.Properties.resize(jsonProperties.size());
      for (size_t i = 0; i < jsonProperties.size(); ++i)
      {
        std::string propertyUUID = jsonProperties[i]["UUID"];
        for (auto &building : baseArea.Buildings)
          if (building.PropertyUUID == propertyUUID)
          {
            baseArea.Properties[i].Buildings.push_back(building);
            break;
          }
        Deserialize(baseArea.Properties[i], json, propertyUUID);
      }
    }

    /// Serialize PrimaryArea
    static void Serialize(const PrimaryArea &primaryArea,
                          nlohmann::json &json,
                          bool isTopObject = true)
    {
      nlohmann::json jsonPrimArea;
      jsonPrimArea["area_id"] = primaryArea.AreaID;
      jsonPrimArea["name"] = primaryArea.Name;
      jsonPrimArea["parent_id"] = primaryArea.DistrictAreaID;
      for (const auto &baseArea : primaryArea.BaseAreas)
      {
        jsonPrimArea["contains"].push_back(baseArea.AreaID);
        Serialize(baseArea, json, false);
      }
      SerializeArea(primaryArea, json, isTopObject, jsonPrimArea,
                    "PrimaryArea");
    }

    /// Deserialize PrimaryArea.
    static void Deserialize(PrimaryArea &primaryArea,
                            const nlohmann::json &json,
                            std::string areaID)
    {
      /// Get corresponding JSON object
      bool isPrimArea = json["Type"] == "PrimaryArea";
      nlohmann::json jsonPrimArea =
          isPrimArea ? json
                     : GetObjectByAttribute("area_id", std::move(areaID),
                                            json["PrimaryAreas"]);
      CheckType("PrimaryArea", jsonPrimArea);

      primaryArea.DistrictAreaID = jsonPrimArea["parent_id"];
      DeserializeArea(primaryArea, jsonPrimArea);
      primaryArea.Name = jsonPrimArea["name"];
      nlohmann::json jsonBaseAreas = jsonPrimArea["contains"];
      for (auto &jsonBaseArea : jsonBaseAreas)
      {
        BaseArea baseArea;
        std::string baseAreaID = jsonBaseArea.get<std::string>();
        Deserialize(baseArea, json, baseAreaID);
        primaryArea.BaseAreas.push_back(baseArea);
      }
    }

    /// Deserialize District.
    static void
    Deserialize(District &district, nlohmann::json &json, std::string areaID)
    {
      /// Get corresponding JSON object
      bool isDistrict = json["Type"] == "District";
      nlohmann::json jsonDistrict =
          isDistrict ? json
                     : GetObjectByAttribute("area_id", std::move(areaID),
                                            json["Districts"]);
      CheckType("District", jsonDistrict);

      DeserializeArea(district, jsonDistrict);
      district.Name = jsonDistrict["name"];
      nlohmann::json jsonPrimAreas = jsonDistrict["contains"];
      for (auto &jsonPrimArea : jsonPrimAreas)
      {
        PrimaryArea primaryArea;
        Deserialize(primaryArea, json, jsonPrimArea.get<std::string>());
        district.PrimaryAreas.push_back(primaryArea);
      }
    }

    /// Serialize District
    static void Serialize(const District &district,
                          nlohmann::json &json,
                          bool isTopObject = true)
    {
      nlohmann::json jsonDistrict;
      jsonDistrict["area_id"] = district.AreaID;
      jsonDistrict["name"] = district.Name;
      for (const auto &primArea : district.PrimaryAreas)
      {
        jsonDistrict["contains"].push_back(primArea.AreaID);
        Serialize(primArea, json, false);
      }
      SerializeArea(district, json, isTopObject, jsonDistrict, "District");
    }

  private:
    /// Serialize (part of) subdivision
    template <typename T>
    static void SerializeArea(const T &t,
                              nlohmann::json &json,
                              bool isTopObject,
                              nlohmann::json &jsonArea,
                              std::string typeName)
    {
      jsonArea["Type"] = typeName;
      SerializeFootprint(t.Footprint, jsonArea);
      if (!isTopObject)
        json[typeName == "Property" ? "Properties" : (typeName + "s")]
            .push_back(jsonArea);
      else
        for (nlohmann::json::iterator it = jsonArea.begin();
             it != jsonArea.end(); ++it)
          json[it.key()] = it.value();
    }

    /// Serialize a footprint.
    static void SerializeFootprint(const Polygon &footprint,
                                   nlohmann::json &json)
    {
      auto jsonFootprint = nlohmann::json::array();
      for (const auto &point : footprint.Vertices)
      {
        nlohmann::json vertex;
        vertex["x"] = point.x;
        vertex["y"] = point.y;
        jsonFootprint.push_back(vertex);
      }
      json["Footprint"] = jsonFootprint;
    }

    /// Deserialize (part of) District, PrimaryArea or BaseArea
    template <typename T>
    static void DeserializeArea(T &subdivision, const nlohmann::json &jsonArea)
    {
      subdivision.AreaID = jsonArea["area_id"];
      Polygon footprint;
      DeserializeFootprint(footprint, jsonArea, "Footprint");
      subdivision.Footprint = footprint;
    }

    /// Deserialize a footprint.
    static void DeserializeFootprint(Polygon &footprint,
                                     const nlohmann::json &json,
                                     const std::string &key = "footprint")
    {
      for (const auto &jsonVertex : json[key])
      {
        Point2D vertex;
        vertex.x = jsonVertex["x"];
        vertex.y = jsonVertex["y"];
        footprint.Vertices.push_back(vertex);
      }
    }

    static nlohmann::json getJsonRoadProps(
        const std::unordered_map<std::string, std::vector<std::string>>
            &properties)
    {
      auto jsonProps = nlohmann::json::object();
      for (const auto &prop : properties)
      {
        nlohmann::json j_vec(prop.second);
        jsonProps[prop.first] = j_vec;
      }
      return jsonProps;
    }
  };

  } // namespace DTCC_BUILDER

#endif
