// SHP I/O
// Anders Logg 2019
// Licensed under the MIT License

#ifndef DTCC_SHP_H
#define DTCC_SHP_H

#include <iostream>
#include <json.hpp>
#include <shapefil.h>
#include <vector>
#include "Geometry.h"
#include "Logging.h"
#include "Polygon.h"

using namespace nlohmann;

namespace DTCC
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
  static void Read(std::vector<Polygon> &polygons,
                   const std::string &fileName,
                   basic_json<> *attributes)
  {
    Info("SHP: Reading polygons from file " + fileName);
    // Open file(s)
    SHPHandle handle = SHPOpen(fileName.c_str(), "r");

    // Get info
    int numEntities, shapeType;
    SHPGetInfo(handle, &numEntities, &shapeType, NULL, NULL);
    Info("SHP: " + str(numEntities) + " entities");
    switch (shapeType)
    {
    case SHPT_POINT:
    case SHPT_POINTZ:
    case SHPT_POINTM:
      Info("SHP: point type");
      break;
    case SHPT_ARC:
    case SHPT_ARCZ:
    case SHPT_ARCM:
      Info("SHP: arc type");
      break;
    case SHPT_POLYGON:
    case SHPT_POLYGONZ:
    case SHPT_POLYGONM:
      Info("SHP: polygon type");
      break;
    case SHPT_MULTIPOINT:
    case SHPT_MULTIPOINTZ:
    case SHPT_MULTIPOINTM:
      Info("SHP: multipoint type");
      break;
    default:
      Info("SHP: unknown type");
    }

    // Check that we have polygon or arc type.
    // TODO: Include SHPT_ARCZ & SHPT_ARCM as well?
    if (shapeType != SHPT_POLYGON && shapeType != SHPT_POLYGONZ &&
        shapeType != SHPT_POLYGONM && shapeType != SHPT_ARC)
      throw std::runtime_error("Shapefile not of relevant type.");

    // Read attributes, if present
    if (attributes != nullptr)
      readAttributes(fileName, numEntities, attributes);

    ReadPolygons(polygons, handle, numEntities);
  }

private:
  /// Remove characters from string.
  /// \param str String to cleanse
  /// \param chars String where all characters to remove should be put
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
      Info("SHP: DBF attributes encoded as ISO-8859-1, converting to UTF-8...");
    else if (codePage != "UTF8")
      Info("SHP: Unknown or unrecognized encoding of DBF attributes, "
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
      attributes->push_back(shapeAttr);
    }
  }

  /// Returns DBF handle.
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

  /// Reads and stores polygons/edges.
  /// \param polygons Vector for polygon storage
  /// \param handle SHP handle
  /// \param numEntities Number of polygons
  static void
  ReadPolygons(std::vector<Polygon> &polygons, SHPInfo *handle, int numEntities)
  {
    for (int i = 0; i < numEntities; i++)
    {
      // Read vertices
      // Get object
      SHPObject *object = SHPReadObject(handle, i);

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

        // Add polygon
        polygons.push_back(polygon);
      }
      else
      {
        // For donut polygons only get the outer hull
        // For multipatch polygons only get the first polygon
        // TODO: handle donut and multipatch polygons correctly

        Polygon polygon;
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
          if (Geometry::PolygonOrientation2D(polygon) == 1)
          {
            polygons.push_back(polygon);
          }
        }
      }
    }
  }
};

} // namespace DTCC

#endif
