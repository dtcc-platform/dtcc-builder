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

/// City represents a collection of buildings.
class City : public Printable
{
public:
  /// name of the city
  std::string name;

  /// Array of buildings
  std::vector<Building> buildings;

  Vector2D origin = Vector2D(0, 0);

  /// Create empty city
  City() = default;
  virtual ~City() {} // make the destructor virtual

  /// Set new origin (subtract offset)
  void set_origin(const Vector2D &origin)
  {
    info("City: Setting new origin to " + str(origin));
    for (auto &building : buildings)
      building.set_origin(origin);
    this->origin = origin;
  }

  /// build search tree (bounding box tree), required for search queries.
  ///
  /// @param rebuild Force rebuild of existing tree if set
  /// @param margin Margin to use for bounding boxes around buildings
  void build_search_tree(bool rebuild = false, double margin = 0.0) const
  {
    // Skip if already built or force rebuild
    if (!bbtree.empty() && !rebuild)
    {
      info("Search tree already built; set rebuild flag to force rebuild.");
      return;
    }

    // Create 2D bounding boxes for all building footprints
    std::vector<BoundingBox2D> bboxes;
    for (const auto &building : buildings)
    {
      BoundingBox2D bbox(building.footprint.vertices, margin);
      bboxes.push_back(bbox);
    }

    // build bounding box tree
    bbtree.build(bboxes);
    debug(str(bbtree));
  }

  // find building containing point (inside footprint), returning -1
  // if the point is not inside any building.
  int find_building(const Vector2D &p) const
  {
    // Check that search tree has been created
    if (bbtree.empty())
    {
      warning("Warning: empty search tree; call build_search_tree()");
      return -1;
    }

    // find candidate buildings from search tree
    std::vector<size_t> indices = bbtree.find(p);

    // Check candidate buildings
    for (const auto index : indices)
    {
      if (Geometry::polygon_contains_2d(buildings[index].footprint, p))
        return index;
    }

    // Point not inside a building
    return -1;
  }

  // 3D version, only using x and y
  int find_building(const Vector3D &p) const
  {
    Vector2D _p(p.x, p.y);
    return find_building(_p);
  }

  // Compute center of city
  Vector2D center() const
  {
    Vector2D c{};
    size_t num_points = 0;
    for (auto const &building : buildings)
    {
      for (auto const &p : building.footprint.vertices)
      {
        c += Vector2D(p);
        num_points += 1;
      }
    }
    c /= static_cast<double>(num_points);
    return c;
  }

  // Compute radius of city (relative to center)
  double radius(const Vector2D &center) const
  {
    double r2_max = 0.0;
    for (auto const &building : buildings)
    {
      for (auto const &p : building.footprint.vertices)
      {
        const double r2 = Geometry::squared_distance_2d(p, center);
        if (r2 > r2_max)
          r2_max = r2;
      }
    }
    return std::sqrt(r2_max);
  }

  /// Pretty-print
  std::string __str__() const override
  {
    return "City " + name + " with " + str(buildings.size()) + " buildings";
  }

private:
  friend class CityBuilder;

  // Bounding box tree used for search queries
  mutable BoundingBoxTree2D bbtree;
};

} // namespace DTCC_BUILDER

#endif
