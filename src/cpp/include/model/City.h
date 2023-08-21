// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_CITY_H
#define DTCC_CITY_H

#include <iomanip>
#include <string>
#include <vector>

#include "BoundingBoxTree.h"
#include "Building.h"
#include "Geometry.h"
#include "Logging.h"
#include "Vector.h"

namespace DTCC_BUILDER
{

/// City represents a collection of Buildings.
class City : public Printable
{
public:
  /// Name of the city
  std::string Name;

  /// Array of buildings
  std::vector<Building> Buildings;

  Point2D Origin = Point2D(0, 0);

  /// Create empty city
  City() = default;
  virtual ~City() {} // make the destructor virtual

  /// Set new origin (subtract offset)
  void SetOrigin(const Point2D &origin)
  {
    info("City: Setting new origin to " + str(origin));
    for (auto &building : Buildings)
      building.SetOrigin(origin);
    Origin = origin;
  }

  /// Build search tree (bounding box tree), required for search queries.
  ///
  /// @param rebuild Force rebuild of existing tree if set
  /// @param margin Margin to use for bounding boxes around buildings
  void BuildSearchTree(bool rebuild = false, double margin = 0.0) const
  {
    // Skip if already built or force rebuild
    if (!bbtree.Empty() && !rebuild)
    {
      info("Search tree already built; set rebuild flag to force rebuild.");
      return;
    }

    // Create 2D bounding boxes for all building footprints
    std::vector<BoundingBox2D> bboxes;
    for (const auto &building : Buildings)
    {
      BoundingBox2D bbox(building.Footprint.Vertices, margin);
      bboxes.push_back(bbox);
    }

    // Build bounding box tree
    bbtree.Build(bboxes);
    debug(str(bbtree));
  }

  // Find building containing point (inside footprint), returning -1
  // if the point is not inside any building.
  int FindBuilding(const Vector2D &p) const
  {
    // Check that search tree has been created
    if (bbtree.Empty())
    {
      warning("Warning: Empty search tree; call BuildSearchTree()");
      return -1;
    }

    // Find candidate buildings from search tree
    std::vector<size_t> indices = bbtree.Find(p);

    // Check candidate buildings
    for (const auto index : indices)
    {
      if (Geometry::PolygonContains2D(Buildings[index].Footprint, p))
        return index;
    }

    // Point not inside a building
    return -1;
  }

  // 3D version, only using x and y
  int FindBuilding(const Vector3D &p) const
  {
    Vector2D _p(p.x, p.y);
    return FindBuilding(_p);
  }

  // Compute center of city
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

  // Compute radius of city (relative to center)
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
  std::string __str__() const override
  {
    return "City " + Name + " with " + str(Buildings.size()) + " buildings";
  }

private:
  friend class CityBuilder;

  // Bounding box tree used for search queries
  mutable BoundingBoxTree2D bbtree;
};

} // namespace DTCC_BUILDER

#endif
