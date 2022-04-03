// Copyright (C) 2020 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_POINT_CLOUD_PROCESSOR_H
#define DTCC_POINT_CLOUD_PROCESSOR_H

#include "KDTreeVectorOfVectorsAdaptor.h"
#include "nanoflann.hpp"

#include "Color.h"
#include "GeoRaster.h"
#include "PointCloud.h"
#include "Timer.h"

namespace DTCC
{

class PointCloudProcessor
{
public:
  // Georaster is assumed to have at least 3 bands where 1,2,3 are R,G,B
  // and in range 0-255
  static void ColorFromImage(PointCloud &pointCloud, const GeoRaster &raster)
  {
    double r, g, b;
    pointCloud.Colors.clear();

    for (auto p : pointCloud.Points)
    {
      Point2D p2d = Point2D(p.x, p.y);
      try
      {
        r = raster(p2d, 1);
        g = raster(p2d, 2);
        b = raster(p2d, 3);
        r /= 255;
        g /= 255;
        b /= 255;
      }
      catch (const std::runtime_error &error)
      {
        // point outside of raster
        r = 0.0;
        g = 0.0;
        b = 0.0;
      }

      pointCloud.Colors.push_back(Color(r, g, b));
    }
  }

  static PointCloud
  ClassificationFilter(const PointCloud &pointCloud,
                       const std::vector<int> &classifications)
  {
    PointCloud outCloud;
    size_t pointCount = pointCloud.Points.size();
    bool has_color = (pointCount == pointCloud.Colors.size());

    if (pointCount != pointCloud.Classifications.size())
    {
      throw std::runtime_error("Point cloud isn't classified");
    }

    for (size_t i = 0; i < pointCount; i++)
    {
      for (auto c : classifications)
      {
        if (pointCloud.Classifications[i] == c)
        {
          outCloud.Points.push_back(pointCloud.Points[i]);
          if (has_color)
          {
            outCloud.Colors.push_back(pointCloud.Colors[i]);
          }
          outCloud.Classifications.push_back(pointCloud.Classifications[i]);
          break;
        }
      }
    }

    // Set bounding box
    outCloud.BoundingBox = pointCloud.BoundingBox;

    return outCloud;
  }

  /// Remove outliers from point cloud by removing all points more than a
  /// given number of standard deviations from the mean for the z-coordinate,
  ///
  /// @param pointCloud The point cloud
  /// @param outlierMargin Number of standard deviations
  //
  // Developer note: This function calls RemoveOutliers() below but needs
  // to do a bit more since also the colors and classifications are affected.
  static void RemoveOutliers(PointCloud &pointCloud, double outlierMargin)
  {
    Info("PointCloudProcessor: Removing outliers...");
    Timer timer("RemoveOutliers");

    // Write heights to file for debbuging
    const bool debugOutliers = true;
    if (debugOutliers)
    {
      std::ofstream f;
      f.open("heights_before.csv");
      for (size_t i = 0; i < pointCloud.Points.size(); i++)
      {
        f << pointCloud.Points[i].z;
        if (i + 1 < pointCloud.Points.size())
          f << ",";
      }
      f.close();
    }

    // Remove outliers from points
    std::vector<size_t> outliers =
        RemoveOutliers(pointCloud.Points, outlierMargin, true);

    // Initialize new colors and classifications
    std::vector<Color> newColors{};
    std::vector<uint8_t> newClassifications{};

    // Copy colors and classifications for all non-outliers
    assert(pointCloud.Colors.size() == pointCloud.Classifications.size());
    size_t k = 0;
    for (size_t i = 0; i < pointCloud.Colors.size(); i++)
    {
      if (k >= outliers.size() || i != outliers[k])
      {
        newColors.push_back(pointCloud.Colors[i]);
        newClassifications.push_back(pointCloud.Classifications[i]);
      }
      else
      {
        k++;
      }
    }

    // Assign new to old
    pointCloud.Colors = newColors;
    pointCloud.Classifications = newClassifications;

    // Write heights to file for debbuging
    if (debugOutliers)
    {
      std::ofstream f;
      f.open("heights_after.csv");
      for (size_t i = 0; i < pointCloud.Points.size(); i++)
      {
        f << pointCloud.Points[i].z;
        if (i + 1 < pointCloud.Points.size())
          f << ",";
      }
      f.close();
    }

    Info("PointCloudProcessor: " + str(outliers.size()) +
         " outliers removed from point cloud");
  }

  /// Remove outliers from vector of points by removing all points more than a
  /// given number of standard deviations from the mean for the z-coordinate.
  ///
  /// @param points The vector of points
  /// @param outlierMargin Number of standard deviations
  /// @return Vector of indices for removed points
  static std::vector<size_t> RemoveOutliers(std::vector<Point3D> &points,
                                            double outlierMargin,
                                            bool verbose = false)
  {
    // Check that we have enough points
    if (points.size() < 3)
      return std::vector<size_t>();

    // Compute min and max
    double min{std::numeric_limits<int>::max()};
    for (const auto p : points)
      min = std::min(min, p.z);
    double max{std::numeric_limits<int>::min()};
    for (const auto p : points)
      max = std::max(max, p.z);

    // Compute mean
    double mean{0};
    for (const auto p : points)
      mean += p.z;
    mean /= points.size();

    // Compute standard deviation
    double std{0};
    for (const auto p : points)
      std += (p.z - mean) * (p.z - mean);
    std /= points.size() - 1;
    std = std::sqrt(std);

    if (verbose)
    {
      Info("PointCloudProcessor: min height = " + str(min) +
           " m (before filtering)");
      Info("PointCloudProcessor: max height = " + str(max) +
           " m (before filtering)");
      Info("PointCloudProcessor: mean height = " + str(mean) + " m");
      Info("PointCloudProcessor: standard deviation = " + str(std) + " m");
    }

    // Remove outliers (can perhaps be implemented more efficiently)
    std::vector<Point3D> newPoints;
    std::vector<size_t> outliers;
    for (size_t i = 0; i < points.size(); i++)
    {
      const Point3D &p = points[i];
      if (std::abs(p.z - mean) <= outlierMargin * std)
      {
        newPoints.push_back(p);
      }
      else
      {
        outliers.push_back(i);
      }
    }
    points = newPoints;

    // Recompute min and max
    min = std::numeric_limits<int>::max();
    for (const auto p : points)
      min = std::min(min, p.z);
    max = std::numeric_limits<int>::min();
    for (const auto p : points)
      max = std::max(max, p.z);

    if (verbose)
    {
      Info("PointCloudProcessor: min height = " + str(min) +
           " m (after filtering)");
      Info("PointCloudProcessor: max height = " + str(max) +
           " m (after filtering)");
    }

    return outliers;
  }
  static std::vector<std::vector<double>>
  KNNNearestNeighbours(std::vector<Point3D> &points, size_t neighbours)
  {
    size_t pc_size = points.size();
    std::vector<std::vector<double>> neighbourDist(pc_size);

    if (neighbours <= 0 or neighbours > pc_size)
    {
      neighbours = pc_size;
    }
    neighbours++; // N neighbours other than ourselves

    typedef KDTreeVectorOfVectorsAdaptor<std::vector<Point3D>, double,
                                         3 /* dims */>
        my_kd_tree_t;
    my_kd_tree_t pc_index(3 /*dim*/, points, 10 /* max leaf */);
    pc_index.index->buildIndex();
    std::vector<size_t> ret_indexes(neighbours);
    std::vector<double> out_dists_sqr(neighbours);
    nanoflann::KNNResultSet<double> resultSet(neighbours);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);

    size_t idx = 0;
    for (auto const &pt : points)
    {
      std::vector<double> query_pt{pt.x, pt.y, pt.z};
      pc_index.query(&query_pt[0], neighbours, &ret_indexes[0],
                     &out_dists_sqr[0]);
      for (size_t i = 1; i < neighbours;
           i++) // start from 1 since 0 is the query point
      {
        neighbourDist[idx].push_back(std::sqrt(out_dists_sqr[i]));
      }
      idx++;
    }

    return neighbourDist;
  }

  /// Finds outliers from vector of points by removing all points more than a
  /// given number of standard deviations from the mean distance to their N
  /// nearest neighbours
  ///
  /// @param points The vector of points
  /// @param neighbours Number of neighbours to consider. If less than 1 or
  /// greater than the number of points in the point cloud use all points
  /// @param outlierMargin Number of standard deviations
  /// @return Vector of indices of outlier points
  static std::vector<size_t>
  StatisticalOutlierFinder(std::vector<Point3D> &points,
                           size_t neighbours,
                           double outlierMargin,
                           bool verbose = false)
  {
    Timer("StatisticalOurlierFinder");
    // Check that we have enough points
    if (points.size() <= neighbours)
      return std::vector<size_t>();

    std::vector<size_t> outliers;

    auto neighbourDist = KNNNearestNeighbours(points, neighbours);
    std::vector<double> u_dist_i;

    for (size_t i = 0; i < points.size(); i++)
    {
      double dsum = 0;
      for (auto &d : neighbourDist[i])
      {
        dsum += d;
      }
      u_dist_i.push_back(dsum / neighbours);
    }

    // Compute mean
    double mean{0};
    for (auto p : u_dist_i)
      mean += p;
    mean /= u_dist_i.size();

    // Compute standard deviation
    double std{0};
    for (auto p : u_dist_i)
      std += (p - mean) * (p - mean);
    std /= u_dist_i.size() - 1;
    std = std::sqrt(std);

    double T = mean + outlierMargin * std;

    // Info("T: " + str(T));
    for (size_t i = 0; i < u_dist_i.size(); i++)
    {
      if (u_dist_i[i] > T)
        outliers.push_back(i);
    }

    return outliers;
  }

  /// Remove outliers from Vector<Point3d> using Statistical Outlier algorithm
  ///
  /// @param points vector of points to filter
  /// @param neighbours Number of neighbours to consider. If less than 1 or
  /// greater than the number of points in the point cloud use all points
  /// @param outlierMargin Number of standard deviations
  /// @param verbose give verbose detail
  static void StatisticalOutlierRemover(std::vector<Point3D> &points,
                                        size_t neighbours,
                                        double outlierMargin,
                                        bool verbose = false)
  {
    Timer("StatisticalOurlierRemover");
    std::vector<size_t> outliers =
        StatisticalOutlierFinder(points, neighbours, outlierMargin, verbose);
    std::vector<Point3D> newPoints;
    size_t k = 0;
    for (size_t i = 0; i < points.size(); i++)
    {
      if (k >= outliers.size() || i != outliers[k])
      {
        newPoints.push_back(points[i]);
      }
      else
      {
        k++;
      }
    }
    points = newPoints;
  }

  /// Remove outliers from PointCloud using Statistical Outlier algorithm
  ///
  /// @param pointCloud The point cloud to filter
  /// @param neighbours Number of neighbours to consider. If less than 1 or
  /// greater than the number of points in the point cloud use all points
  /// @param outlierMargin Number of standard deviations
  /// @param verbose give verbose detail
  static void StatisticalOutlierRemover(PointCloud &pointCloud,
                                        size_t neighbours,
                                        double outlierMargin,
                                        bool verbose = false)
  {
    // points to remove
    std::vector<size_t> outliers = StatisticalOutlierFinder(
        pointCloud.Points, neighbours, outlierMargin, verbose);

    std::vector<Point3D> newPoints;
    std::vector<Color> newColors{};
    std::vector<uint8_t> newClassifications{};
    // Copy points, colors and classifications for all non-outliers
    // assert(pointCloud.Colors.size() == pointCloud.Classifications.size());
    bool hasColor = false;
    bool hasClass = false;
    if (pointCloud.Colors.size() > 0)
      hasColor = true;
    if (pointCloud.Classifications.size() > 0)
      hasClass = true;
    size_t k = 0;
    for (size_t i = 0; i < pointCloud.Points.size(); i++)
    {
      if (k >= outliers.size() || i != outliers[k])
      {
        newPoints.push_back(pointCloud.Points[i]);

        if (hasColor)
          newColors.push_back(pointCloud.Colors[i]);
        if (hasClass)
          newClassifications.push_back(pointCloud.Classifications[i]);
      }
      else
      {
        k++;
      }
    }
    // Assign new to old
    pointCloud.Points = newPoints;

    if (hasColor)
      pointCloud.Colors = newColors;
    if (hasClass)
      pointCloud.Classifications = newClassifications;
  }
};

} // namespace DTCC

#endif
