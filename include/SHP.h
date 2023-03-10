// SHP I/O
// Copyright (C) 2019, 2020 Anders Logg, Anton J Olsson
// Anders Logg, Vasilis Naserentin 2019
// Licensed under the MIT License

#ifndef DTCC_SHP_H
#define DTCC_SHP_H

#include "Geometry.h"
#include "Logging.h"
#include "Polygon.h"
#include "Utils.h"
#include <iostream>
#include <json.hpp>
#include <shapefil.h>
#include <vector>

using namespace nlohmann;

namespace DTCC_BUILDER
{

class SHP
{
public:
  /// Read polygons/edges and (possibly) attributes from SHP file. Note that the
  /// corresponding .shx and .dbf files must also be present in the same
  /// directory.
  ///
  /// \param polygons The vector to put the polygons in
  /// \param fileName The SHP filename
  /// \param attributes A JSON object to put possible attributes in
  /// \param handleDuplicateUUIDs Append unique number to duplicate UUIDs
  static void Read(std::vector<Polygon> &polygons,
                   const std::string &fileName,
                   std::vector<std::string> *UUIDs = nullptr,
                   std::vector<int> *entityID = nullptr,
                   json *attributes = nullptr,
                   bool handleDuplicateUUIDs = true)
  {
    info("SHP: Reading polygons from file " + fileName);
    // Open file(s)
    SHPHandle handle = SHPOpen(fileName.c_str(), "r");
    DBFHandle handleD = DBFOpen(fileName.c_str(), "rb");
    // DBFHandle handleD = getDBFHandle(fileName);

    if (UUIDs != nullptr && entityID == nullptr)
    {
      error("entityID can't be null if reading UUIDs");
    }

    // Get info
    int numEntities, shapeType;
    SHPGetInfo(handle, &numEntities, &shapeType, nullptr, nullptr);
    info("SHP: " + str(numEntities) + " entities");
    switch (shapeType)
    {
    case SHPT_POINT:
    case SHPT_POINTZ:
    case SHPT_POINTM:
      info("SHP: point type");
      break;
    case SHPT_ARC:
    case SHPT_ARCZ:
    case SHPT_ARCM:
      info("SHP: arc type");
      break;
    case SHPT_POLYGON:
    case SHPT_POLYGONZ:
    case SHPT_POLYGONM:
      info("SHP: polygon type");
      break;
    case SHPT_MULTIPOINT:
    case SHPT_MULTIPOINTZ:
    case SHPT_MULTIPOINTM:
      info("SHP: multipoint type");
      break;
    default:
      info("SHP: unknown type");
    }

    // Check that we have polygon or arc type.
    // TODO: Include SHPT_ARCZ & SHPT_ARCM as well?
    if (shapeType != SHPT_POLYGON && shapeType != SHPT_POLYGONZ &&
        shapeType != SHPT_POLYGONM && shapeType != SHPT_ARC)
      throw std::runtime_error("Shapefile not of relevant type.");

    // Read attributes, if present
    if (attributes != nullptr)
      readAttributes(fileName, numEntities, attributes);

    ReadPolygons(polygons, handle, numEntities, shapeType, attributes, handleD,
                 UUIDs, entityID, handleDuplicateUUIDs);
  }

private:
  /// Remove characters from string.
  ///
  /// \param str String to cleanse
  /// \param chars String with all characters to remove
  /// \return The cleansed string
  static std::string removeChars(std::string &str, const std::string &chars)
  {
    str.erase(std::remove_if(str.begin(), str.end(),
                             [chars](unsigned char x) {
                               return chars.find(x) != std::string::npos;
                             }),
              str.end());
    return str;
  }

  /// Get all field/attribute names.
  ///
  /// \param handle The DBF handle.
  /// \param numFields Number of fields in the DBF file
  /// \param codePage The encoding of the DBF file
  /// \return Vector of field names
  static std::vector<std::string>
  getFieldNames(DBFHandle handle, int numFields, const std::string &codePage)
  {
    std::vector<std::string> fieldNames;
    for (int i = 0; i < numFields; ++i)
    {
      char fieldName[12];
      DBFGetFieldInfo(handle, i, fieldName, nullptr, nullptr);
      std::string fieldNameStr(fieldName);
      if (codePage == "ISO88591")
        fieldNameStr = Utils::Iso88591ToUtf8(fieldNameStr);
      fieldNames.emplace_back(fieldNameStr);
    }
    return fieldNames;
  }

  /// Read and store all attributes/fields.
  ///
  /// \param fileName The SHP filename
  /// \param numEntities Number of polygons
  /// \param attributes JSON object to store the attributes in
  static void readAttributes(const std::string &fileName,
                             int numEntities,
                             basic_json<> *attributes)
  {
    DBFHandle handle = getDBFHandle(fileName);
    std::string codePage;
    if (handle->pszCodePage != nullptr)
      codePage = DBFGetCodePage(handle);
    codePage = removeChars(codePage, " -");
    std::transform(codePage.begin(), codePage.end(), codePage.begin(),
                   ::toupper);
    if (codePage == "ISO88591")
      info("SHP: DBF attributes encoded as ISO-8859-1, converting to "
           "UTF-8...");
    else if (codePage != "UTF8")
      info("SHP: Unknown or unrecognized encoding of DBF attributes, "
           "characters may "
           "not be "
           "displayed correctly");

    int numFields = DBFGetFieldCount(handle);
    std::vector<std::string> fieldNames =
        getFieldNames(handle, numFields, codePage);

    for (int i = 0; i < numEntities; ++i)
    {
      json shapeAttr = json({});
      for (int j = 0; j < numFields; ++j)
      {
        std::string attribute = DBFReadStringAttribute(handle, i, j);
        if (codePage == "ISO88591")
          attribute = Utils::Iso88591ToUtf8(attribute);
        shapeAttr[fieldNames[j]] = attribute;
      }
      (*attributes)["attributes"].push_back(shapeAttr);
    }
  }

  /// Return DBF handle.
  ///
  /// \param fileName SHP filename.
  /// \return The DBF handle.
  static DBFHandle getDBFHandle(const std::string &fileName)
  {
    DBFHandle dbfHandle;
    std::string dbfName = fileName;
    dbfName.replace(fileName.length() - 3, 3, "dbf");
    dbfHandle = DBFOpen(dbfName.c_str(), "rb");
    return dbfHandle;
  }

  /// Append -N to UUID in case of duplicates.
  /// \param UUIDs Vector with all UUIDs fetched so far
  /// \param UUID The UUID to compare with other UUIDs
  /// \return Unique UUID
  static std::string MakeUUIDUnique(const std::vector<std::string> *UUIDs,
                                    const std::string &UUID)
  {
    int duplicates = 0;
    for (const std::string &oldUUID : *UUIDs)
    {
      if (oldUUID.find(UUID) != std::string::npos)
        duplicates++;
    }
    return UUID + (duplicates > 0 ? "-" + std::to_string(duplicates) : "");
  }

  /// Add an entity's data to corresponding containers.
  /// \param polygons Polygon vector
  /// \param UUIDs UUID vector
  /// \param entityID 1-based indices vector
  /// \param index Index of entity to be added
  /// \param UUID UUID of entity
  /// \param polygon Entity's polygon
  /// \param handleDuplicateUUIDs Append unique number to duplicate UUIDs
  static void AddData(std::vector<Polygon> &polygons,
                      std::vector<std::string> *UUIDs,
                      std::vector<int> *entityID,
                      int index,
                      const std::string &UUID,
                      const Polygon &polygon,
                      bool handleDuplicateUUIDs = true)
  {
    // Add polygon
    polygons.push_back(polygon);
    // Add entityID and UUID
    if (UUIDs != nullptr)
    {
      UUIDs->push_back(handleDuplicateUUIDs ? MakeUUIDUnique(UUIDs, UUID)
                                            : UUID);
      entityID->push_back(index + 1); // 1-based index
    }
  }

  static std::string
  getUUID(DBFHandle handleD, const std::vector<std::string> *UUIDs, int index)
  {
    // Open DFB to read UUID
    std::string test;
    if (UUIDs != nullptr)
    {
      test = DBFReadStringAttribute(handleD, index, 0);
      if (test.empty())
        test = Utils::CreateUUID() + "::DTCCGenerated";
    }
    return test;
  }

  /// Read and store polygons/edges.
  ///
  /// \param polygons Vector for polygon storage
  /// \param handle SHP handle
  /// \param numEntities Number of polygons
  /// \param handleDuplicateUUIDs Append unique number to duplicate UUIDs
  static void ReadPolygons(std::vector<Polygon> &polygons,
                           SHPInfo *handle,
                           int numEntities,
                           int shapeType,
                           json *attributes,
                           DBFHandle handleD,
                           std::vector<std::string> *UUIDs,
                           std::vector<int> *entityID,
                           bool handleDuplicateUUIDs = true)
  {
    for (int i = 0; i < numEntities; i++)
    {
      // Read vertices
      std::string test = getUUID(handleD, UUIDs, i);

      // Get object
      SHPObject *object = SHPReadObject(handle, i);

      if (attributes != nullptr)
      {
        for (int j = 0; j < object->nParts; ++j)
        {
          (*attributes)["polyToAttribute"].push_back(i);
        }
      }

      // Get vertices
      if (object->nParts == 1)
      {
        // Create empty polygon
        Polygon polygon;

        for (int j = 0; j < object->nVertices; j++)
        {
          const double x = object->padfX[j];
          const double y = object->padfY[j];
          Vector2D p(x, y);
          polygon.Vertices.push_back(p);
        }
        AddData(polygons, UUIDs, entityID, i, test, polygon,
                handleDuplicateUUIDs);
      }
      else
      {
        // For donut polygons only get the outer hull
        // For multipatch polygons save each part as its own polygon
        // TODO: handle donut and multipatch polygons correctly

        int start;
        int end;
        for (int part = 0; part < object->nParts; part++)
        {
          Polygon polygon;
          start = object->panPartStart[part];
          if (part + 1 == object->nParts)
          {
            end = object->nVertices;
          }
          else
          {
            end = object->panPartStart[part + 1];
          }

          for (int j = start; j < end; j++)
          {
            const double x = object->padfX[j];
            const double y = object->padfY[j];
            Vector2D p(x, y);
            polygon.Vertices.push_back(p);
          }
          if (Geometry::PolygonOrientation2D(polygon) == 1 ||
              shapeType == SHPT_ARC)
            AddData(polygons, UUIDs, entityID, i, test, polygon,
                    handleDuplicateUUIDs);
        }
      }
    }
  }
};

} // namespace DTCC_BUILDER

#endif
