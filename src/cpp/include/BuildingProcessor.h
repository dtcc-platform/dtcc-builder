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
  static double PointCoverage(const Building &building, double tileSize = 1.0)
  {
    // Estimate what percentage of a building roof is covered by the point cloud
    // If much less than 1 then that indicates that there is a problem with the
    // pointcloud data for that building.
    Timer("BuildingProcessor::PointCoverage");
    auto bbox = BoundingBox2D(building.Footprint);
    std::vector<BoundingBox2D> tiles;
    for (double x = bbox.P.x; x < bbox.Q.x; x += tileSize)
    {
      for (double y = bbox.P.y; y < bbox.Q.y; y += tileSize)
      {
        auto tile =
            BoundingBox2D(Point2D(x, y), Point2D(x + tileSize, y + tileSize));
        Polygon tilePoly;
        tilePoly.Vertices.push_back(Point2D(x, y));
        tilePoly.Vertices.push_back(Point2D(x + tileSize, y));
        tilePoly.Vertices.push_back(Point2D(x + tileSize, y + tileSize));
        tilePoly.Vertices.push_back(Point2D(x, y + tileSize));

        if (Geometry::Intersects2D(building.Footprint, tilePoly))
        {
          tiles.push_back(tile);
        }
      }
    }

    BoundingBoxTree2D tileTree;
    tileTree.Build(tiles);

    std::set<size_t> collisionSet;
    for (const auto &point : building.RoofPoints)
    {
      auto p = Point2D(point.x, point.y);
      auto containingTile = tileTree.Find(p);

      if (containingTile.size() > 0)
      {
        collisionSet.insert(containingTile[0]);
      }
    }

    return static_cast<double>(collisionSet.size()) / tiles.size();
  }
};

} // namespace DTCC_BUILDER

#endif
