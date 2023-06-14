// Copyright (C) 2020-2022 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_POINT_CLOUD_PROCESSOR_H
#define DTCC_POINT_CLOUD_PROCESSOR_H

#include <Eigen/SVD>
#include <fstream>
#include <math.h>
#include <random>

#include "KDTreeVectorOfVectorsAdaptor.h"
#include "nanoflann.hpp"

#include "Timer.h"
#include "model/Color.h"
#include "model/GeoRaster.h"
#include "model/Point.h"
#include "model/PointCloud.h"
#include "model/Vector.h"

template <typename T> int sign(T val) { return (T(0) < val) - (val < T(0)); }

namespace DTCC_BUILDER
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
    info("PointCloudProcessor: Removing outliers...");
    Timer timer("RemoveOutliers");

    // Write heights to file for debbuging
    const bool debugOutliers = false;
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
        FindGlobalOutliers(pointCloud.Points, outlierMargin);

    FilterPointCloud(pointCloud, outliers);

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

    info("PointCloudProcessor: " + str(outliers.size()) +
         " outliers removed from point cloud");
  }

  /// Find index of outlier from vector of points more than a
  /// given number of standard deviations from the mean for the z-coordinate.
  ///
  /// @param points The vector of points
  /// @param outlierMargin Number of standard deviations
  /// @return Vector of indices for removed points
  static std::vector<size_t>
  FindGlobalOutliers(const std::vector<Point3D> &points, double outlierMargin)
  {
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
    std::vector<size_t> outliers;

    for (size_t i = 0; i < points.size(); i++)
    {
      if (std::abs(points[i].z - mean) > outlierMargin * std)
      {
        outliers.push_back(i);
      }
    }

    return outliers;
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
      info("PointCloudProcessor: min height = " + str(min) +
           " m (before filtering)");
      info("PointCloudProcessor: max height = " + str(max) +
           " m (before filtering)");
      info("PointCloudProcessor: mean height = " + str(mean) + " m");
      info("PointCloudProcessor: standard deviation = " + str(std) + " m");
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
      info("PointCloudProcessor: min height = " + str(min) +
           " m (after filtering)");
      info("PointCloudProcessor: max height = " + str(max) +
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
    Timer("StatisticalOurtierFinder");
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

    // info("T: " + str(T));
    for (size_t i = 0; i < u_dist_i.size(); i++)
    {
      if (u_dist_i[i] > T)
        outliers.push_back(i);
    }

    return outliers;
  }

  /// Returns the Distance to the K nearest neighbors for each point in
  /// points
  static std::vector<std::vector<double>>
  KNNNearestNeighboursDist(std::vector<Point3D> &points, size_t neighbours)
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

  /// Returns the index of the K nearest neighbors for each point in
  /// points
  static std::vector<std::vector<size_t>>
  KNNNearestNeighboursIdx(std::vector<Point3D> &points, size_t neighbours)
  {
    size_t pc_size = points.size();
    std::vector<std::vector<size_t>> neighbourIdx(pc_size);

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
        neighbourIdx[idx].push_back(ret_indexes[i]);
      }
      idx++;
    }

    return neighbourIdx;
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
    Timer("StatisticalOurtierRemover");
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

    FilterPointCloud(pointCloud, outliers);
  }

  static void RANSAC_OutlierRemover(PointCloud &pointCloud,
                                    double distanceThreshold,
                                    size_t iterations = 100)
  {
    Timer("RANSAC_OutlierRemover");
    std::vector<size_t> outliers =
        RANSAC_OutlierFinder(pointCloud.Points, distanceThreshold, iterations);
    FilterPointCloud(pointCloud, outliers);
  }

  static void RANSAC_OutlierRemover(std::vector<Point3D> &points,
                                    double distanceThreshold,
                                    size_t iterations = 100)
  {
    Timer("RANSAC_OutlierRemover");

    auto outliers = RANSAC_OutlierFinder(points, distanceThreshold, iterations);
    if (outliers.size() == 0)
      return;
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

  static std::vector<size_t> RANSAC_OutlierFinder(std::vector<Point3D> &points,
                                                  double distanceThreshold,
                                                  size_t iterations = 100)
  {

    std::vector<size_t> outliers;
    if (points.size() < 9)
      return outliers;
    std::vector<size_t> best_outliers(points.size(), 0);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    std::default_random_engine generator(seed);
    std::uniform_int_distribution<size_t> distribution(0, points.size() - 1);
    auto randIdx = std::bind(distribution, generator);
    Point3D pt1, pt2, pt3;
    Vector3D v1, v2, v3;
    size_t idx1, idx2, idx3;
    double k;

    for (size_t i = 0; i < iterations; i++)
    {
      idx1 = randIdx();
      idx2 = randIdx();
      // idx1, 2 and 3 must be different
      while (idx2 == idx1)
      {
        idx2 = randIdx();
      }
      idx3 = randIdx();
      while (idx3 == idx1 || idx3 == idx2)
      {
        idx3 = randIdx();
      }
      pt1 = points[idx1];
      pt2 = points[idx2];
      pt3 = points[idx3];
      v1 = Vector3D(pt2.x - pt1.x, pt2.y - pt1.y, pt2.z - pt1.z);
      v2 = Vector3D(pt3.x - pt1.x, pt3.y - pt1.y, pt3.z - pt1.z);
      v3 = v1.Cross(v2);
      v3 /= v3.Magnitude();
      if (isnan(v3.x)) // all three points are in a line
      {
        continue;
      }

      k = v3.Dot(pt2);
      outliers.clear();
      for (size_t j = 0; j < points.size(); j++)
      {
        auto ptPlaneDist = std::abs(v3.Dot(points[j]) - k) / v3.Magnitude();
        if (ptPlaneDist > distanceThreshold)
        {

          outliers.push_back(j);
        }
      }
      if (outliers.size() < best_outliers.size())
      {
        best_outliers = outliers;
      }
    }
    return best_outliers;
  }

  // Removes selected points from point cloud
  static void FilterPointCloud(PointCloud &pointCloud,
                               std::vector<size_t> ptsToRemove)
  {
    std::sort(ptsToRemove.begin(), ptsToRemove.end());
    std::vector<Point3D> newPoints;
    std::vector<Vector3D> newNormals;
    std::vector<Color> newColors{};
    std::vector<uint8_t> newClassifications{};
    std::vector<uint16_t> newIntensities{};
    std::vector<uint8_t> newScanFlags{};

    bool hasNormals = pointCloud.Normals.size() > 0;
    bool hasColor = pointCloud.Colors.size() > 0;
    bool hasClass = pointCloud.Classifications.size() > 0;
    bool hasIntensity = pointCloud.Intensities.size() > 0;
    bool hasScanflags = pointCloud.ScanFlags.size() > 0;

    size_t k = 0;
    for (size_t i = 0; i < pointCloud.Points.size(); i++)
    {
      if (k >= ptsToRemove.size() || i != ptsToRemove[k])
      {
        newPoints.push_back(pointCloud.Points[i]);

        if (hasNormals)
          newNormals.push_back(pointCloud.Normals[i]);
        if (hasColor)
          newColors.push_back(pointCloud.Colors[i]);
        if (hasClass)
          newClassifications.push_back(pointCloud.Classifications[i]);
        if (hasIntensity)
          newIntensities.push_back(pointCloud.Intensities[i]);
        if (hasScanflags)
          newScanFlags.push_back(pointCloud.ScanFlags[i]);
      }
      else
      {
        k++;
      }
    }
    // Assign new to old
    pointCloud.Points = newPoints;

    if (hasNormals)
      pointCloud.Normals = newNormals;
    if (hasColor)
      pointCloud.Colors = newColors;
    if (hasClass)
      pointCloud.Classifications = newClassifications;
    if (hasIntensity)
      pointCloud.Intensities = newIntensities;
    if (hasScanflags)
      pointCloud.ScanFlags = newScanFlags;
  }

  static std::pair<uint8_t, uint8_t> parseScanFlag(uint8_t flag)
  {
    uint8_t returnNumber = flag & 7;
    uint8_t numReturns = (flag >> 3) & 7;
    return std::pair<uint8_t, uint8_t>(returnNumber, numReturns);
  }

  static uint8_t packScanFlag(uint8_t returnNumber, uint8_t numReturns)
  {
    return (returnNumber & 7) | ((numReturns & 7) << 3);
  }

  static void NaiveVegetationFilter(PointCloud &pointCloud)
  {
    // Remove some points that might be vegetation
    std::vector<size_t> pointsToRemove;
    if (pointCloud.ScanFlags.size() != pointCloud.Points.size())
    {
      warning("Scan flags not set. No vegetation filtering");
      return;
    }
    bool hasClassification = false;
    if (pointCloud.Classifications.size() == pointCloud.Points.size())
      hasClassification = true;

    for (size_t i = 0; i < pointCloud.Points.size(); i++)
    {
      auto scanFlag = parseScanFlag(pointCloud.ScanFlags[i]);
      uint8_t classification = 1;
      if (hasClassification)
        classification = pointCloud.Classifications[i];
      // not last point of several

      if ((classification >= 3 && classification <= 5) || // classified as veg
          (classification < 2 &&                          // not classified
           (scanFlag.first != scanFlag.second))           // not last
      )
      {
        pointsToRemove.push_back(i);
      }
    }
    FilterPointCloud(pointCloud, pointsToRemove);
  }

  static void EstimateNormalsKNN(PointCloud &pointCloud, size_t neighbours)
  {
    auto normals = EstimateNormalsKNN(pointCloud.Points, neighbours);
    pointCloud.Normals = normals;
  };

  static std::vector<Vector3D> EstimateNormalsKNN(std::vector<Point3D> points,
                                                  size_t neighbours)
  {
    std::vector<Vector3D> normals;
    auto neigboursIdx = KNNNearestNeighboursIdx(points, neighbours);
    size_t idx = 0;
    Eigen::RowVector3d dir(0.0, 0.0, 1.0);
    for (auto const &query_pt : points)
    {
      size_t found = neigboursIdx[idx].size();
      if (found < 3) // not enough neighbours to estimate normal
      {
        normals.push_back(Vector3D(0, 0, 0));
        idx++;
        continue;
      }

      Eigen::MatrixXd neighbors(found, 3);
      for (int i = 0; i < found; i++)
      {
        auto pt = points[neigboursIdx[idx][i]];
        neighbors(i, 0) = pt.x - query_pt.x;
        neighbors(i, 1) = pt.y - query_pt.y;
        neighbors(i, 2) = pt.z - query_pt.z;
      }
      Eigen::RowVector3d normal;
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(neighbors, Eigen::ComputeThinV);
      Eigen::MatrixXd V = svd.matrixV();
      for (int l = 0; l < 3; l++)
      {
        normal[l] = V(l, 2);
      }
      normal *= sign(normal.dot(dir));
      auto n = Vector3D(normal[0], normal[1], normal[2]);
      // Info("Normal: " + str(n));
      normals.push_back(n);

      idx++;
    }

    return normals;
  }
};

} // namespace DTCC_BUILDER

#endif
