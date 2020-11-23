// SHP I/O
// Anders Logg 2019

#ifndef DTCC_SHP_H
#define DTCC_SHP_H

#include <iostream>
#include <shapefil.h>
#include <vector>

#include "Polygon.h"
#include "Geometry.h"
#include "Logging.h"

namespace DTCC
{

class SHP
{
public:
  // Read polygons from SHP file. Note that the corresponding
  // .shx and .dbf files must also be present in the same directory.
  static void Read(std::vector<Polygon> &polygons,
                   std::vector<std::string> &UUIDs,
                   std::string fileName)
  {
    Info("SHP: Reading polygons from file " + fileName);
    // Open file(s)
    SHPHandle handle = SHPOpen(fileName.c_str(), "r");
    DBFHandle handleD = DBFOpen(fileName.c_str(), "r");

    // Get info
    int numEntities, shapeType;
    SHPGetInfo(handle, &numEntities, &shapeType, nullptr, nullptr);
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

    // Check that we have polygon type
    if (shapeType != SHPT_POLYGON && shapeType != SHPT_POLYGONZ && shapeType != SHPT_POLYGONM)
      throw std::runtime_error("Shapefile not of polygon type.");

    // Read footprints
    for (int i = 0; i < numEntities; i++)
    {
      // Get UUID
      const char *test = DBFReadStringAttribute(handleD, i, 0);
      UUIDs.push_back(test);
      // std::cout<<UUIDs[i]<<std::endl;

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

        //Polygon polygon;
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
            end =  object->panPartStart[part + 1];
          }

          for (int j = start; j < end; j++)
          {
            const double x = object->padfX[j];
            const double y = object->padfY[j];
            Vector2D p(x, y);
            polygon.Vertices.push_back(p);
          }
          if  (Geometry::PolygonOrientation2D(polygon) == 1)
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
