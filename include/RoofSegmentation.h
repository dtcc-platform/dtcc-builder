// Copyright (C) 2020 Anton J Olsson
// Licensed under the MIT License

#ifndef DTCC_ROOFSEGMENTATION_H
#define DTCC_ROOFSEGMENTATION_H

#include <random>
#include <set>

#include "KDTreeVectorOfVectorsAdaptor.h"
#include "Logging.h"
#include "Point.h"
#include "PointCloudProcessor.h"
#include "Polygon.h"
#include "Utils.h"
#include "Vector.h"
#include "datamodel/Building.h"

namespace DTCC
{
class RoofSegmentation
{
public:
  static void SegmentRoofs(std::vector<Building> &buildings,
                           double max_radius = 2,
                           double normal_angle_threshold = 2,
                           double seed_ratio = 0.1)
  {
    for (auto &building : buildings)
    {
      auto regions = RegionGrowingSegmentation(
          building.RoofPoints, max_radius, normal_angle_threshold, seed_ratio);
    }
  }

  static std::vector<std::vector<size_t>>
  RegionGrowingSegmentation(const std::vector<Point3D> &points,
                            double max_radius = 2,
                            double normal_angle_threshold = 2,
                            double seed_ratio = 0.1)
  {
    std::vector<std::vector<size_t>> regions;

    auto normals = PointCloudProcessor::EstimateNormalsKNN(points, 8);

    typedef KDTreeVectorOfVectorsAdaptor<std::vector<Point3D>, double,
                                         3 /* dims */>
        my_kd_tree_t;
    my_kd_tree_t pc_index(3, points, 10 /* max leaf */);
    max_radius = max_radius * max_radius;
    normal_angle_threshold = normal_angle_threshold * M_PI / 180;

    size_t seed_point_size = (size_t)round(seed_ratio * points.size());
    std::set<size_t> seed_points_idx;
    while (seed_points_idx.size() < seed_point_size)
    {
      seed_points_idx.insert(std::rand() % points.size());
    }

    std::set<size_t> assigned_points;
    for (size_t idx : seed_points_idx)
    {
      std::set<size_t> S;
      S.insert(idx);
      std::vector<size_t> R;
      while (S.size() > 0)
      {
        size_t p_idx = *S.begin();
        S.erase(S.begin());

        auto p = points[p_idx];
        std::vector<double> query_pt{p.x, p.y, p.z};
        auto p_normal = normals[p_idx];
        auto p_neighbors = pc_index.radiusQuery(&query_pt[0], max_radius);
        for (auto const &ind_pt : p_neighbors)
        {
          size_t q_idx = ind_pt.first;
          if (assigned_points.count(q_idx) == 0)
          {
            double angle_between = p_normal.AngleBetween(normals[q_idx]);
            if (angle_between < normal_angle_threshold)
            {
              R.push_back(q_idx);
              S.insert(q_idx);
              assigned_points.insert(q_idx);
            }
          }
        }
      }
      if (R.size() > 0)
      {
        regions.push_back(R);
      }
    }
    return regions;
  }
};
} // namespace DTCC

#endif