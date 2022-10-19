// Copyright (C) 2020 Anders Logg, Vasilis Naserentin and Orfeas Eleftheriou
// Licensed under the MIT License

#ifndef DTCC_CITY_JSON_H
#define DTCC_CITY_JSON_H

//#include "JSON.h"
#include "Point.h"
#include "Logging.h"
#include <utility>
#include <vector>
#include <unordered_map>
#include <iostream>

namespace DTCC_BUILDER
{

/// CityObject is a type that stores parsed data from the CityJSON format
class CityObject : public Printable
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

  /// Various attributes of the CityObject.
  /// Currently not all fields are parsed from the json format.
  /// More info available at the official docs:
  /// https://portal.opengeospatial.org/files/?artifact_id=47842
  struct Attributes
  {
    /// Retrieved measured height
    double MeasuredHeight;

    /// Retrieved roof type
    std::string RoofType;

    /// Stories above ground
    int StoreysAboveGround;

    /// Description of this attribute
    std::string Description;

    /// Name of this attribute
    std::string Name;
    // TODO: Add all the field from the "Attributes" type of CityJSON based on
    // its documentation

    Attributes() : MeasuredHeight{0}, StoreysAboveGround{0} {};
  };

  struct Geometry
  {

    /// Declaring types found inside Geometry field of CityJSON
    /// Currently only Solid, GeometryInstance and MultiSurface are supported.
    enum GeometryType
    {
      Solid,
      GeometryInstance,
      MultiSurface
      // TODO: Add all the types from Geometry type of CityJSON based on its
      // documentation
    };

    /// Boundary type acts as a wrapper for Boundaries ID located inside the Boundaries array in the CityJSON format
    struct Boundary
    {
      /// Contains all IDs for the current Boundary
      std::vector<uint> BoundariesIDs;

      friend std::ostream &operator<<(std::ostream &os,
                                      const Boundary &boundary)
      {
        std::string boundariesStr;
        for(auto it=boundary.BoundariesIDs.begin();it!=boundary.BoundariesIDs.end();++it)
        {
            boundariesStr.append(str(*it) + ",");
        }
        return os << boundariesStr;
      }
    };

    // Declaring members of Geometry
    /// LOD level of this Geometry
    uint LOD;

    /// GeometryType for this Geometry
    GeometryType Type;

    /// All boundaries for this geometry
    std::vector<Boundary> Boundaries;

    Geometry() : LOD{0}, Type{GeometryType::Solid} {};

    //TODO: Move this to source file in next iteration
    /*static std::unordered_map<std::string, GeometryType> const GeometryTypesStringTable = {
              {"Solid", GeometryType::Solid},
              {"GeometryInstance", GeometryType::GeometryInstance},
              {"MultiSurface",GeometryType::MultiSurface}};*/

    static std::string GeometryTypeToString(GeometryType typeToConvert)
    {
      switch (typeToConvert)
      {
      case GeometryType::Solid:
        return "Solid";
      case GeometryType::GeometryInstance:
        return "GeometryInstance";
      case GeometryType::MultiSurface:
        return "MultiSurface";
      default:
        return "Solid";
      }
    }
  };


public:

  /// ID of this Object
  std::string ID;

  /// CityJSON attributes for this Object
  Attributes ObjectAttributes;

  /// CityJSON Object type for this Object
  CityObjectType ObjectType;

  /// Retrieved geometry for this Object
  Geometry ObjectGeometry;

  //TODO: Add geographical extent, children and parents based on CityJSON's docs

  /// Create ampty CityObject
  CityObject() : ID("N/A") {}

  /// Create city object with given ID
  ///
  /// @param newID the id of this object
  explicit CityObject(std::string newId) : ID(std::move(newId)) {};

  /// Pretty-print
  virtual std::string __str__() const
  {
    std::string BuildingStr;

    // boundaries print
    std::string boundariesStr;
    for (uint i = 0; i < ObjectGeometry.Boundaries.size(); i++)
    {
      for (uint j = 0; j < ObjectGeometry.Boundaries[i].BoundariesIDs.size();
           j++)
      {
        boundariesStr.append(
            std::to_string(ObjectGeometry.Boundaries[i].BoundariesIDs[j]) +
            ","); // TODO: refactor this
      }
      boundariesStr.erase(boundariesStr.size() - 1); // remove last comma
      boundariesStr.append("|");
    }
    // Remove last character from the string since it contains garbage data
    // coming from the for loop before this
    boundariesStr.erase(boundariesStr.size() - 1);

    BuildingStr.append("ID:" + ID + " - Boundaries:");
    BuildingStr.append(boundariesStr);
    return BuildingStr;
  }

  friend std::ostream &operator<<(std::ostream &os, const CityObject &obj)
  {
    return os << obj.__str__() << std::endl;
  }

};

/// CityJSON acts as a data container that stores all the relevant information from CityJSON's files
class CityJSON
{

public:

  /// Contains all the vertices for the parsed geometry
  std::vector<Point3D> Vertices;

  /// Contains all the city objects of CityJSON file
  std::vector<CityObject> CityObjects;

  /// Type of the CityJSON file
  std::string Type;

  /// Version of the CityJSON file
  std::string Version;

public:

    /// Create Empty CityJSON
    CityJSON() : Type("CityJSON"), Version("1.0") {}
};

} // namespace DTCC_BUILDER
#endif
