// Copyright (C) 2022 Dag Wästberg
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

typedef std::vector<std::vector<size_t>> Segments;
class RoofSegmentation
{
public:
  static void SegmentRoofs(std::vector<Building> &buildings,
                           double max_radius = 2,
                           double normal_angle_threshold = 2,
                           size_t min_points = 10,
                           double seed_ratio = 0.1)
  {
    for (auto &building : buildings)
    {
      auto regions = RegionGrowingSegmentation(building.RoofPoints, max_radius,
                                               normal_angle_threshold,
                                               min_points, seed_ratio);
      building.RoofSegments = regions;
    }
  }

  static void SegmentRoofsSearchOptimal(std::vector<Building> &buildings,
                                        double max_radius = 2,
                                        double seed_ratio = 0.1)
  {
    for (auto &building : buildings)
    {
      double best_angle = SearchOptimateRegionSegmentationAngle(
          building.RoofPoints, max_radius, seed_ratio);
      auto regions = RegionGrowingSegmentation(building.RoofPoints, max_radius,
                                               best_angle, seed_ratio);
      building.RoofSegments = regions;
    }
  }

  static double
  SearchOptimateRegionSegmentationAngle(const std::vector<Point3D> &points,
                                        double max_radius = 2,
                                        size_t min_points = 10,
                                        double seed_ratio = 0.1)
  {
    double best_score = -1;
    double prev_score = -1;
    double best_angle = -1;
    double step = 0.25;
    double start = 0.5;
    double end = 10;
    double a = start;
    while (a < end)
    {
      auto roof_segments = RegionGrowingSegmentation(points, max_radius, a,
                                                     min_points, seed_ratio);
      double score = scoreRegion(roof_segments, points);
      if (score > best_score)
      {
        best_score = score;
        best_angle = a;
      }
      if (score > prev_score)
        a += step;
      else
        a *= 2;
      prev_score = score;
    }
    return best_angle;
  }

  static Segments RegionGrowingSegmentation(const std::vector<Point3D> &points,
                                            double max_radius = 2,
                                            double normal_angle_threshold = 2,
                                            size_t min_points = 12,
                                            double seed_ratio = 0.1)
  {
    Segments regions;

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
        // Info("S: " + str(S.size()));
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
            // Info("angle_between: " + str(angle_between));
            if (angle_between < normal_angle_threshold)
            {
              R.push_back(q_idx);
              S.insert(q_idx);
              assigned_points.insert(q_idx);
            }
          }
        }
      }
      if (R.size() >= min_points)
      {
        regions.push_back(R);
      }
    }
    return regions;
  }

private:
  static double scoreRegion(Segments regions,
                            const std::vector<Point3D> &points)
  {
    double score;
    double segmented_points = 0;
    for (auto const &region : regions)
    {
      segmented_points += region.size();
    }
    double segment_ratio = segmented_points / points.size();
    score = regions.size() * regions.size() * segment_ratio;
    return score;
  }
};
} // namespace DTCC

#endif