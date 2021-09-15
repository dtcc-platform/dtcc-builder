// Copyright (C) 2020 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_POINT_CLOUD_PROCESSOR_H
#define DTCC_POINT_CLOUD_PROCESSOR_H

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
      double r,g,b;
      pointCloud.Colors.clear();

      for (auto p: pointCloud.Points) {
        Point2D p2d = Point2D(p.x,p.y);
        try
        {
          r = raster(p2d,1);
          g = raster(p2d,2);
          b = raster(p2d,3);
          r /= 255;
          g /= 255;
          b /= 255;
        } catch (const std::runtime_error& error)
        {
          // point outside of raster
          r = 0.0;
          g = 0.0;
          b = 0.0;
        }

        pointCloud.Colors.push_back(Color(r,g,b));
      }

  }

  static PointCloud ClassificationFilter(const PointCloud &pointCloud, const std::vector<int> &classifications)
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
};

} // namespace DTCC

#endif
