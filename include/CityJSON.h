// Representation of a CityJSON
// Copyright (C) 2020 Anders Logg, Vasilis Naserentin, Orfeas Eleftheriou.

#ifndef DTCC_CITYJSON_H
#define DTCC_CITYJSON_H

#include "JSON.h"
#include "Point.h"
#include <vector>
#include <unordered_map>
#include <iostream>
/*
 {
  "type": "CityJSON",
  "version": "1.0",
  "CityObjects": {},
  "vertices": []
}

{
 "CityObjects": {
         "attributes": {},
         "geometry": [
           { "boundaries": [], "lod": 1, "semantics":{"surfaces":[],
"values":[]}, "type": "Solid"
             }
           ],
           "type": "Building"
 },
  "type": "CityJSON",
  "version": "1.0",
  "vertices": [
  ]
}

*/
namespace DTCC
{
// struct Attributes {
//    double MeasuredHeight;
//};
// struct Semantics {
//    std::vector<nlohmann::json> Surfaces;
//    std::vector<nlohmann::json> Values;
//};
// struct _Geometry {
//    std::vector<nlohmann::json> Boundaries;
//    int64_t Lod;
//    Semantics _Semantics;
//    std::string Type;
//};
// struct CityObjects {
//    Attributes _Attributes;
//    std::vector<_Geometry> Geometry;
//    std::string Type;
//};

class CityObject
{
public:
  enum CityObjectType
  {
    Building,
    Transportation,
    TINRelief,
    WaterBody
    // TODO: Add all the City Object types from CityJson's documentation
  };

  struct Attributes
  {
    double MeasuredHeight;
    std::string RoofType;
    int StoreysAboveGround;
    std::string Description;
    std::string Name;
    // TODO: Add all the field from the "Attributes" type of CityJSON based on
    // its documentation
  };

  struct Geometry
  {

    // Declaring types found inside Geometry field of CityJSON
    enum GeometryType
    {
      Solid,
      GeometryInstance,
      MultiSurface
      // TODO: Add all the types from Geometry type of CityJSON based on its
      // documentation
    };

    struct Boundary
    {
      std::vector<uint> BoundariesIDs;

      friend std::ostream& operator<<(std::ostream& os, const Boundary& boundary) 
      {
        std::string boundariesStr;
        for(auto it=boundary.BoundariesIDs.begin();it!=boundary.BoundariesIDs.end();++it)
        {
            boundariesStr.append(*it +",");
        }
        return os << boundariesStr;
      }
    };

    // Declaring members of Geometry
    uint LOD;
    GeometryType Type;
    std::vector<Boundary> Boundaries;

    //TODO: Move this to source file in next iteration
    /*static std::unordered_map<std::string, GeometryType> const GeometryTypesStringTable = {
              {"Solid", GeometryType::Solid},
              {"GeometryInstance", GeometryType::GeometryInstance},
              {"MultiSurface",GeometryType::MultiSurface}};*/


  };


private:
  std::string ID;
  Attributes attributes;
  CityObjectType objectType;
  Geometry objectGeometry;
  // TODO: Add geographical extent, children and parents based on CityJSON's
  // docs

public:

  CityObject() : ID("N/A") {}
  CityObject(std::string newId) : ID(newId) {};

  //Public getters

  inline std::string GetID() { return ID; }

  //Public setters 
  inline void SetObjectType(CityObjectType NewObjectType) { objectType = NewObjectType;}
  inline void SetAttributes(Attributes NewAttributes) { attributes=NewAttributes; }
  inline void SetObjectGeometry(Geometry NewGeometry) { objectGeometry = NewGeometry; }

  friend std::ostream &operator<<(std::ostream &os, const CityObject &obj) 
  {
    //boundaries print
    std::string boundariesStr;
    for(uint i=0;i<obj.objectGeometry.Boundaries.size();i++)
    {
      for(uint j=0;j<obj.objectGeometry.Boundaries[i].BoundariesIDs.size();j++)
      {
        boundariesStr.append(std::to_string(obj.objectGeometry.Boundaries[i].BoundariesIDs[j]) +","); // TODO: refactor this
      }
      boundariesStr.erase(boundariesStr.size()-1); //remove last comma 
      boundariesStr.append("|");
    }
    //Remove last character from the string since it contains garbage data coming from the for loop before this
    boundariesStr.erase(boundariesStr.size()-1);
    return os << "ID:" << obj.ID << " - Boundaries: " << boundariesStr/*obj.objectGeometry.Boundaries*/ << std::endl;
  }

};

class CityJSON
{
  // public:
  //
  //    std::vector<std::string> UUID;
  //    CityObjects _CityObjects;
  //    std::string Type;
  //    std::string Version;
  //    nlohmann::json Vertices;
  //    // Create empty CityJSON
  //    CityJSON() : Type("CityJSON"), Version("1.0"), Vertices{} {}
  //
  //    void printInfo()
  //    {
  //        Info("Information of CityJSON object");
  //        Info(Type);
  //        Info(Version);
  //        Info("Included UUIDs: ");
  //        for (size_t i = 0; i < UUID.size(); i++)
  //        {
  //            Info(UUID[i]);
  //        }
  //    }

private:

  std::vector<Point3D> vertices;
  std::vector<CityObject> cityObjects;
  std::string type;
  std::string version;

  void storeVertices(nlohmann::json jsonFile) 
  {
    uint verticesNum = jsonFile["vertices"].size();
    for(uint i=0;i<verticesNum;i++)
    {
      //debug
      std::cout << jsonFile["vertices"][i]<<std::endl;
      auto vertex = jsonFile["vertices"][i];

      assert(vertex.size()==3 && "Attempted to read non 3D vertex from json file");

      vertices.push_back(Point3D(jsonFile["vertices"][i][0],
                                 jsonFile["vertices"][i][1],
                                 jsonFile["vertices"][i][2]));

    }
  }

  CityObject::Attributes getAttributesFromJson(nlohmann::json jsonFile, const std::string& cityObjectID)
  {
    CityObject::Attributes attributes=CityObject::Attributes();
    std::cout << jsonFile["CityObjects"][cityObjectID]["attributes"]<<std::endl;
    attributes.MeasuredHeight = jsonFile["CityObjects"][cityObjectID]["attributes"]["measuredHeight"];
    //TODO: Get more fields for attributes here
    
    return attributes;
  }

  CityObject::Geometry getGeometryFromJson(nlohmann::json jsonFile, const std::string& cityObjectID)
  {
    CityObject::Geometry geometry = CityObject::Geometry();
    auto jsonGeometryArray = jsonFile["CityObjects"][cityObjectID]["geometry"][0];
    std::cout<<jsonGeometryArray<<std::endl;
    
    //Storing lod settings
    geometry.LOD = jsonGeometryArray["lod"];
    //std::cout << "LOD:"<<geometry.LOD<<std::endl;

    //Storing boundaries
    uint boundariesArraySize = jsonGeometryArray["boundaries"][0].size();
    for(uint i=0;i<boundariesArraySize;i++)
    {
      auto boundariesJson = jsonGeometryArray["boundaries"][0][i];
      for(auto boundaryIter = boundariesJson.begin();boundaryIter!=boundariesJson.end();++boundaryIter)
      {
        //std::cout << "boundary iter:" << *boundaryIter << std::endl;
        auto internalJsonArray = boundariesJson[0];

        CityObject::Geometry::Boundary newBoundary;

        for (auto it = internalJsonArray.begin(); it != internalJsonArray.end();
             ++it)
        {
          //std::cout << *it << std::endl;
          newBoundary.BoundariesIDs.push_back(*it);
        }

        geometry.Boundaries.push_back(newBoundary);
      }
    }
    
    //Storing building type
    std::string tempGeomType = jsonGeometryArray["type"];
    //TODO: correlate this to string-enum
    if(tempGeomType.compare("Solid"))
    {
        geometry.Type=CityObject::Geometry::Solid;
    }
    else
    {
    //TODO: this should be replaced by built-in functionality in struct
    }
    

    return geometry;
  }

  void storeCityObjects(nlohmann::json jsonFile) 
  {
    auto cityObjectsField = jsonFile["CityObjects"].get<nlohmann::json::object_t>();
    for(auto cityObjectIterator=cityObjectsField.begin();cityObjectIterator!=cityObjectsField.end();++cityObjectIterator)
    {
        std::cout<<cityObjectIterator->first<<std::endl; //getting City Object ID here

        //Creating a new object in stack, we're going to store the information here and later store this into cityObjects vector
        CityObject NewObj = CityObject(cityObjectIterator->first);
        
        auto jsonCityObj =
            jsonFile["CityObjects"][cityObjectIterator->first];

        //Getting attributes
        CityObject::Attributes attributes = getAttributesFromJson(jsonFile,NewObj.GetID());
        NewObj.SetAttributes(attributes);

        //Getting geometry
        CityObject::Geometry geometry = getGeometryFromJson(jsonFile, NewObj.GetID());
        NewObj.SetObjectGeometry(geometry);

        //Getting type
        // TODO: move this whole impl to cityobj class function
        std::string tempStr = jsonCityObj["type"];
        if(tempStr.compare("Building"))
        {
            NewObj.SetObjectType(CityObject::Building);
        }

        cityObjects.push_back(NewObj);
    }
  }

public:

    CityJSON() : type("CityJSON"), version("1.0") {}
    
    /**
     * Tries to parse a json file of format CityJSON and to fill out the corresponding fields.
     * TODO: Convert this from a constructor to a separate function and move the whole implementation in JSON.h
     */
    CityJSON(nlohmann::json jsonFile) : CityJSON()
    {

        storeVertices(jsonFile);
        storeCityObjects(jsonFile);
     /* std::vector<aiVector3D> vertices;
      uint num = jsonFile["vertices"].size();
      for (uint i = 0; i < num; i++)
      {
        std::cout << jsonFile["vertices"][i] << std::endl;
        auto vertex = jsonFile["vertices"][i];
        std::vector<double> tempVertexArray;
        for (auto it = vertex.begin(); it != vertex.end(); ++it)
        {
          tempVertexArray.push_back((double)*it);
        }
        vertices.push_back(aiVector3D(tempVertexArray[0], tempVertexArray[1],
                                      tempVertexArray[2]));
      }*/


    }

    inline std::vector<CityObject> getCityObjects() const { return cityObjects; }

};

} // namespace DTCC
#endif
