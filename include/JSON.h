// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_JSON_H
#define DTCC_JSON_H

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

#include "Parameters.h"
#include "BoundingBox.h"
#include "Grid.h"
#include "Mesh.h"
#include "Surface.h"
#include "GridField.h"
#include "GridVectorField.h"
#include "CityModel.h"
#include "CityJSON.h"
#include "Color.h"
#include "ColorMap.h"
#include "Road.h"

namespace DTCC
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
      parameters.DataDirectory = ToString("DataDirectory", json);
      parameters.X0 = ToDouble("X0", json);
      parameters.Y0 = ToDouble("Y0", json);
      parameters.XMin = ToDouble("XMin", json);
      parameters.YMin = ToDouble("YMin", json);
      parameters.XMax = ToDouble("XMax", json);
      parameters.YMax = ToDouble("YMax", json);
      parameters.AutoDomain =ToBool("AutoDomain", json);
      parameters.DomainHeight = ToDouble("DomainHeight", json);
      parameters.ElevationModelResolution =
          ToDouble("ElevationModelResolution", json);
      parameters.MeshResolution = ToDouble("MeshResolution", json);
      parameters.MinimalBuildingDistance = ToDouble("MinimalBuildingDistance", json);
      parameters.MinimalVertexDistance =
          ToDouble("MinimalVertexDistance", json);
      parameters.FlatGround = ToBool("FlatGround", json);
      parameters.GroundSmoothing = ToInt("GroundSmoothing", json);
      parameters.Debug = ToBool("Debug", json);
    };

    /// Serialize Parameters
    static void Serialize(const Parameters& parameters, nlohmann::json& json)
    {
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
      json["ElevationModelResolution"] = parameters.ElevationModelResolution;
      json["MeshResolution"] = parameters.MeshResolution;
      json["MinimalBuildingDistance"] = parameters.MinimalBuildingDistance;
      json["MinimalVertexDistance"] = parameters.MinimalVertexDistance;
      json["FlatGround"] = parameters.FlatGround;
      json["GroundSmoothing"] = parameters.GroundSmoothing;
      json["Debug"] = parameters.Debug;
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

    /// Serialize Mesh2D
    static void Serialize(const Mesh2D& mesh, nlohmann::json& json)
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

    /// Serialize Mesh3D
    static void Serialize(const Mesh3D& mesh, nlohmann::json& json)
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
    }

    /// Deserialize Surface3D
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
      const auto jsonCells = json["Cells"];
      surface.Cells.resize(jsonCells.size() / 2);
      for (size_t i = 0; i < surface.Cells.size(); i++)
      {
        surface.Cells[i].v0 = jsonCells[2*i];
        surface.Cells[i].v1 = jsonCells[2*i + 1];
      }
    }

    /// Serialize Surface3D
    static void Serialize(const Surface2D& surface, nlohmann::json& json)
    {
      auto jsonVertices = nlohmann::json::array();
      for (const auto p: surface.Vertices)
      {
        jsonVertices.push_back(p.x);
        jsonVertices.push_back(p.y);
      }
      auto jsonCells = nlohmann::json::array();
      for (const auto T: surface.Cells)
      {
        jsonCells.push_back(T.v0);
        jsonCells.push_back(T.v1);
      }
      json["Type"] = "Surface2D";
      json["Vertices"] = jsonVertices;
      json["Cells"] = jsonCells;
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
      const auto jsonCells = json["Cells"];
      surface.Cells.resize(jsonCells.size() / 3);
      for (size_t i = 0; i < surface.Cells.size(); i++)
      {
        surface.Cells[i].v0 = jsonCells[3*i];
        surface.Cells[i].v1 = jsonCells[3*i + 1];
        surface.Cells[i].v2 = jsonCells[3*i + 2];
      }
    }

    /// Serialize Surface3D
    static void Serialize(const Surface3D& surface, nlohmann::json& json)
    {
      auto jsonVertices = nlohmann::json::array();
      for (const auto p: surface.Vertices)
      {
        jsonVertices.push_back(p.x);
        jsonVertices.push_back(p.y);
        jsonVertices.push_back(p.z);
      }
      auto jsonCells = nlohmann::json::array();
      for (const auto T: surface.Cells)
      {
        jsonCells.push_back(T.v0);
        jsonCells.push_back(T.v1);
        jsonCells.push_back(T.v2);
      }
      json["Type"] = "Surface3D";
      json["Vertices"] = jsonVertices;
      json["Cells"] = jsonCells;
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

    /// Serialize GridField2D
    static void Serialize(const GridField2D& field, nlohmann::json& json)
    {
      auto jsonGrid = nlohmann::json::object();
      Serialize(field.Grid, jsonGrid);
      json["Type"] = "GridField2D";
      json["Grid"] = jsonGrid;
      json["Values"] = field.Values;
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

    /// Serialize GridField3D
    static void Serialize(const GridField3D& field, nlohmann::json& json)
    {
      auto jsonGrid = nlohmann::json::object();
      Serialize(field.Grid, jsonGrid);
      json["Type"] = "GridField3D";
      json["Grid"] = jsonGrid;
      json["Values"] = field.Values;
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

    /// Serialize GridVectorField2D
    static void Serialize(const GridVectorField2D& field, nlohmann::json& json)
    {
      auto jsonGrid = nlohmann::json::object();
      Serialize(field.Grid, jsonGrid);
      json["Type"] = "GridVectorField2D";
      json["Grid"] = jsonGrid;
      json["Values"] = field.Values;
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

    /// Serialize GridVectorField3D
    static void Serialize(const GridVectorField3D& field, nlohmann::json& json)
    {
      auto jsonGrid = nlohmann::json::object();
      Serialize(field.Grid, jsonGrid);
      json["Type"] = "GridVectorField3D";
      json["Grid"] = jsonGrid;
      json["Values"] = field.Values;
    }

    // FIXME: Add separate serialize/deserialize for Building and use here.

    /// Serialize CityModel
    static void Deserialize(CityModel &cityModel, const nlohmann::json& json)
    {
      CheckType("CityModel", json);
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
        cityModel.Buildings[i].GroundHeight = jsonBuilding["GroundHeight"];
       cityModel.Buildings[i].UUID = jsonBuilding["UUID"];
      }
    }

    /// Deserialize CityModel
    static void Serialize(const CityModel &cityModel, nlohmann::json& json)
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
        jsonBuilding["Height"] = building.Height;
        jsonBuilding["GroundHeight"] = building.GroundHeight;
        jsonBuilding["UUID"] = building.UUID;
        jsonBuilding["debugD"] = building.debugID;
        jsonBuildings.push_back(jsonBuilding);
      }
      json["Type"] = "CityModel";
      json["Buildings"] = jsonBuildings;
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

    /// Deserialize Road
    static void Deserialize(Road &road, const nlohmann::json& json)
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
      {
        road.Edges[i] = jsonEdges[i];
      }
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
    }

    /// Serialize Road
    static void Serialize(const Road &road, nlohmann::json& json)
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
      for (const auto e: road.Edges)
      {
        jsonEdges.push_back(e);
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
      // Add to RoadNetwork object
      json["RoadNetwork"] = jsonRoadNetwork;
      // Add Type parameter
      json["Type"] = "RoadNetwork";
    }
  };

} // namespace DTCC

#endif
