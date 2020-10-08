// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_CITY_MODEL_H
#define DTCC_CITY_MODEL_H

#include <iomanip>
#include <string>
#include <vector>

#include "Building.h"
#include "Geometry.h"
#include "Vector.h"
#include "BoundingBoxTree.h"
#include "Logging.h"

namespace DTCC
{

  class CityModel : public Printable
  {
  public:
    // List of buildings
    std::vector<Building> Buildings;

    // Build search tree
    void BuildSearchTree() const
    {
      // Create bounding boxes for all building footprints
      std::vector<BoundingBox2D> bboxes;
      for (const auto& building: Buildings)
      {
        BoundingBox2D bbox(building.Footprint);
        bboxes.push_back(bbox);
      }

      // Build bounding box tree
      bbtree.Build(bboxes);
      std::cout << bbtree << std::endl;
      Progress(str(bbtree));
    }

    // Find building containing point (inside footprint), returning -1
    // if the point is not inside any building.
    int FindBuilding(const Vector2D &p) const
    {
      // Check that search tree has been created
      if (bbtree.Nodes.empty())
      {
        Warning("Warning: Missing search tree; call BuildSearchTree()");
        return -1;
      }

      // Find candidate buildings from search tree
      std::vector<size_t> indices = bbtree.Find(p);

      // Check candidate buildings
      for (const auto index: indices)
      {
        if (Geometry::PolygonContains2D(Buildings[index].Footprint, p))
          return index;
      }

      // Point not inside a building
      return -1;
    }

    // Compute center of city model
    Point2D Center() const
    {
      Vector2D c{};
      size_t numPoints = 0;
      for (auto const &building : Buildings)
      {
        for (auto const &p : building.Footprint.Vertices)
        {
          c += Vector2D(p);
          numPoints += 1;
        }
      }
      c /= static_cast<double>(numPoints);
      return c;
    }

    // Compute radius of city model (relative to center)
    double Radius(const Vector2D &center) const
    {
      double r2max = 0.0;
      for (auto const &building : Buildings)
      {
        for (auto const &p : building.Footprint.Vertices)
        {
          const double r2 = Geometry::SquaredDistance2D(p, center);
          if (r2 > r2max)
            r2max = r2;
        }
      }
      return std::sqrt(r2max);
    }

    /// Pretty-print
    std::string __str__() const
    {
      return "CityModel with " + str(Buildings.size()) + " buildings";
    }

  private:

    // Bounding box tree (used by find building)
    mutable BoundingBoxTree2D bbtree;

  };

} // namespace DTCC

#endif
