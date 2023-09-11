// Copyright (C) 2022 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_BUILDING_PROCESSOR_H
#define DTCC_BUILDING_PROCESSOR_H

#include <set>

#include "BoundingBox.h"
#include "BoundingBoxTree.h"
#include "Geometry.h"
#include "Timer.h"
#include "model/Building.h"
#include "model/Polygon.h"

namespace DTCC_BUILDER
{

class BuildingProcessor
{
public:
  static double point_coverage(const Building &building, double tile_size = 1.0)
  {
    // Estimate what percentage of a building roof is covered by the point cloud
    // If much less than 1 then that indicates that there is a problem with the
    // pointcloud data for that building.
    Timer("BuildingProcessor::point_coverage");
    auto bbox = BoundingBox2D(building.footprint);
    std::vector<BoundingBox2D> tiles;
    for (double x = bbox.P.x; x < bbox.Q.x; x += tile_size)
    {
      for (double y = bbox.P.y; y < bbox.Q.y; y += tile_size)
      {
        auto tile = BoundingBox2D(Vector2D(x, y),
                                  Vector2D(x + tile_size, y + tile_size));
        Polygon tile_poly;
        tile_poly.vertices.push_back(Vector2D(x, y));
        tile_poly.vertices.push_back(Vector2D(x + tile_size, y));
        tile_poly.vertices.push_back(Vector2D(x + tile_size, y + tile_size));
        tile_poly.vertices.push_back(Vector2D(x, y + tile_size));

        if (Geometry::intersects_2d(building.footprint, tile_poly))
        {
          tiles.push_back(tile);
        }
      }
    }

    BoundingBoxTree2D tile_tree;
    tile_tree.build(tiles);

    std::set<size_t> collision_set;
    for (const auto &point : building.roof_points)
    {
      auto p = Vector2D(point.x, point.y);
      auto containingTile = tile_tree.find(p);

      if (containingTile.size() > 0)
      {
        collision_set.insert(containingTile[0]);
      }
    }

    return static_cast<double>(collision_set.size()) / tiles.size();
  }
};

} // namespace DTCC_BUILDER

#endif
