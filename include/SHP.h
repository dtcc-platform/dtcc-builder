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
  // Read polygons and (possibly) attributes from SHP file. Note that the
  // corresponding .shx and .dbf files must also be present in the same
  // directory.
  static void Read(std::vector<Polygon> &polygons,
                   const std::string &fileName,
                   basic_json<> *attributes)
  {
    Info("SHP: Reading polygons from file " + fileName);
    // Open file(s)
    SHPHandle handle = SHPOpen(fileName.c_str(), "r");

    DBFHandle dbfHandle = nullptr;
    if (attributes != nullptr)
      dbfHandle = getDBFHandle(fileName);

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

    ReadPolygons(polygons, handle, numEntities, dbfHandle, attributes);
  }

private:
  static DBFHandle getDBFHandle(const std::string &fileName)
  {
    DBFHandle dbfHandle;
    std::string dbfName = fileName;
    dbfName.replace(fileName.length() - 3, 3, "dbf");
    dbfHandle = DBFOpen(dbfName.c_str(), "rb");
    return dbfHandle;
  }

  static void
  ReadAttributes(DBFHandle handle, basic_json<> *attributes, int shapeIndex)
  {
    static int codeFieldIndex, categoryFieldIndex = -1;
    if (codeFieldIndex < 0)
      codeFieldIndex = DBFGetFieldIndex(handle, "KKOD");
    int code = DBFReadIntegerAttribute(handle, shapeIndex, codeFieldIndex);
    if (categoryFieldIndex < 0)
      categoryFieldIndex = DBFGetFieldIndex(handle, "KATEGORI");
    std::string category = std::string(
        DBFReadStringAttribute(handle, shapeIndex, categoryFieldIndex));
    json shapeAttr = json({});
    shapeAttr["KKOD"] = code;
    shapeAttr["KATEGORI"] = category;
    attributes->push_back(shapeAttr);
  }

  static void ReadPolygons(std::vector<Polygon> &polygons,
                           SHPInfo *handle,
                           int numEntities,
                           DBFHandle dbfHandle,
                           basic_json<> *attributes)
  {
    for (int i = 0; i < numEntities; i++)
    {
      if (dbfHandle != nullptr && attributes != nullptr)
      {
        ReadAttributes(dbfHandle, attributes, i);
      }

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
