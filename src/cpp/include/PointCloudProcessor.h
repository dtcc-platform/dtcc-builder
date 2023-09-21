// Copyright (C) 2020-2022 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_POINT_CLOUD_PROCESSOR_H
#define DTCC_POINT_CLOUD_PROCESSOR_H

#include <Eigen/SVD>
#include <fstream>
#include <iso646.h>
#include <math.h>
#include <random>

#include "KDTreeVectorOfVectorsAdaptor.h"
#include "nanoflann.hpp"

#include "Timer.h"
#include "model/Color.h"
#include "model/GeoRaster.h"
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
  static void color_from_image(PointCloud &point_cloud, const GeoRaster &raster)
  {
    double r, g, b;
    point_cloud.colors.clear();

    for (auto p : point_cloud.points)
    {
      Vector2D p2d = Vector2D(p.x, p.y);
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

      point_cloud.colors.push_back(Color(r, g, b));
    }
  }

  static PointCloud
  classification_filter(const PointCloud &point_cloud,
                        const std::vector<int> &classifications)
  {
    PointCloud out_cloud;
    size_t point_count = point_cloud.points.size();
    bool has_color = (point_count == point_cloud.colors.size());

    if (point_count != point_cloud.classifications.size())
    {
      throw std::runtime_error("Point cloud isn't classified");
    }

    for (size_t i = 0; i < point_count; i++)
    {
      for (auto c : classifications)
      {
        if (point_cloud.classifications[i] == c)
        {
          out_cloud.points.push_back(point_cloud.points[i]);
          if (has_color)
          {
            out_cloud.colors.push_back(point_cloud.colors[i]);
          }
          out_cloud.classifications.push_back(point_cloud.classifications[i]);
          break;
        }
      }
    }

    // Set bounding box
    out_cloud.bounding_box = point_cloud.bounding_box;

    return out_cloud;
  }

  /// Remove outliers from point cloud by removing all points more than a
  /// given number of standard deviations from the mean for the z-coordinate,
  ///
  /// @param point_cloud The point cloud
  /// @param outlier_margin Number of standard deviations
  //
  // Developer note: This function calls remove_outliers() below but needs
  // to do a bit more since also the colors and classifications are affected.
  static void remove_outliers(PointCloud &point_cloud, double outlier_margin)
  {
    info("Removing outliers...");
    Timer timer("remove_outliers");

    // Write heights to file for debbuging
    const bool debug_outliers = false;
    if (debug_outliers)
    {
      std::ofstream f;
      f.open("heights_before.csv");
      for (size_t i = 0; i < point_cloud.points.size(); i++)
      {
        f << point_cloud.points[i].z;
        if (i + 1 < point_cloud.points.size())
          f << ",";
      }
      f.close();
    }

    // Remove outliers from points
    std::vector<size_t> outliers =
        find_global_outliers(point_cloud.points, outlier_margin);

    filter_point_cloud(point_cloud, outliers);

    // Write heights to file for debbuging
    if (debug_outliers)
    {
      std::ofstream f;
      f.open("heights_after.csv");
      for (size_t i = 0; i < point_cloud.points.size(); i++)
      {
        f << point_cloud.points[i].z;
        if (i + 1 < point_cloud.points.size())
          f << ",";
      }
      f.close();
    }

    info(str(outliers.size()) + " outliers removed from point cloud");
  }

  /// find index of outlier from vector of points more than a
  /// given number of standard deviations from the mean for the z-coordinate.
  ///
  /// @param points The vector of points
  /// @param outlier_margin Number of standard deviations
  /// @return Vector of indices for removed points
  static std::vector<size_t>
  find_global_outliers(const std::vector<Vector3D> &points,
                       double outlier_margin)
  {
    // Compute mean
    double mean{0};
    for (const auto &p : points)
      mean += p.z;
    mean /= points.size();

    // Compute standard deviation
    double std{0};
    for (const auto &p : points)
      std += (p.z - mean) * (p.z - mean);
    std /= points.size() - 1;
    std = std::sqrt(std);
    std::vector<size_t> outliers;

    for (size_t i = 0; i < points.size(); i++)
    {
      if (std::abs(points[i].z - mean) > outlier_margin * std)
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
  /// @param outlier_margin Number of standard deviations
  /// @return Vector of indices for removed points
  static std::vector<size_t> remove_outliers(std::vector<Vector3D> &points,
                                             double outlier_margin,
                                             bool verbose = false)
  {
    // Check that we have enough points
    if (points.size() < 3)
      return std::vector<size_t>();

    // Compute min and max
    double min{std::numeric_limits<int>::max()};
    for (const auto &p : points)
      min = std::min(min, p.z);
    double max{std::numeric_limits<int>::min()};
    for (const auto &p : points)
      max = std::max(max, p.z);

    // Compute mean
    double mean{0};
    for (const auto &p : points)
      mean += p.z;
    mean /= points.size();

    // Compute standard deviation
    double std{0};
    for (const auto &p : points)
      std += (p.z - mean) * (p.z - mean);
    std /= points.size() - 1;
    std = std::sqrt(std);

    if (verbose)
    {
      info("min height = " + str(min) + " m (before filtering)");
      info("max height = " + str(max) + " m (before filtering)");
      info("mean height = " + str(mean) + " m");
      info("Standard deviation = " + str(std) + " m");
    }

    // Remove outliers (can perhaps be implemented more efficiently)
    std::vector<Vector3D> new_points;
    std::vector<size_t> outliers;
    for (size_t i = 0; i < points.size(); i++)
    {
      const Vector3D &p = points[i];
      if (std::abs(p.z - mean) <= outlier_margin * std)
      {
        new_points.push_back(p);
      }
      else
      {
        outliers.push_back(i);
      }
    }
    points = new_points;

    // Recompute min and max
    min = std::numeric_limits<int>::max();
    for (const auto &p : points)
      min = std::min(min, p.z);
    max = std::numeric_limits<int>::min();
    for (const auto &p : points)
      max = std::max(max, p.z);

    if (verbose)
    {
      info("min height = " + str(min) + " m (after filtering)");
      info("max height = " + str(max) + " m (after filtering)");
    }

    return outliers;
  }
  static std::vector<std::vector<double>>
  knn_nearest_neighbours(std::vector<Vector3D> &points, size_t neighbours)
  {
    size_t pc_size = points.size();
    std::vector<std::vector<double>> neighbour_dist(pc_size);

    if (neighbours <= 0 or neighbours > pc_size)
    {
      neighbours = pc_size;
    }
    neighbours++; // N neighbours other than ourselves

    typedef KDTreeVectorOfVectorsAdaptor<std::vector<Vector3D>, double,
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
        neighbour_dist[idx].push_back(std::sqrt(out_dists_sqr[i]));
      }
      idx++;
    }

    return neighbour_dist;
  }

  /// Finds outliers from vector of points by removing all points more than a
  /// given number of standard deviations from the mean distance to their N
  /// nearest neighbours
  ///
  /// @param points The vector of points
  /// @param neighbours Number of neighbours to consider. If less than 1 or
  /// greater than the number of points in the point cloud use all points
  /// @param outlier_margin Number of standard deviations
  /// @return Vector of indices of outlier points
  static std::vector<size_t>
  statistical_outlier_finder(std::vector<Vector3D> &points,
                             size_t neighbours,
                             double outlier_margin,
                             bool verbose = false)
  {
    Timer("StatisticalOurtierFinder");
    // Check that we have enough points
    if (points.size() <= neighbours)
      return std::vector<size_t>();

    std::vector<size_t> outliers;

    auto neighbour_dist = knn_nearest_neighbours(points, neighbours);
    std::vector<double> u_dist_i;

    for (size_t i = 0; i < points.size(); i++)
    {
      double dsum = 0;
      for (auto &d : neighbour_dist[i])
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

    double T = mean + outlier_margin * std;

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
  knn_nearest_neighbours_dist(std::vector<Vector3D> &points, size_t neighbours)
  {
    size_t pc_size = points.size();
    std::vector<std::vector<double>> neighbour_dist(pc_size);

    if (neighbours <= 0 or neighbours > pc_size)
    {
      neighbours = pc_size;
    }
    neighbours++; // N neighbours other than ourselves

    typedef KDTreeVectorOfVectorsAdaptor<std::vector<Vector3D>, double,
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
        neighbour_dist[idx].push_back(std::sqrt(out_dists_sqr[i]));
      }
      idx++;
    }

    return neighbour_dist;
  }

  /// Returns the index of the K nearest neighbors for each point in
  /// points
  static std::vector<std::vector<size_t>>
  knn_nearest_neighbours_idx(std::vector<Vector3D> &points, size_t neighbours)
  {
    size_t pc_size = points.size();
    std::vector<std::vector<size_t>> neighbour_idx(pc_size);

    if (neighbours <= 0 or neighbours > pc_size)
    {
      neighbours = pc_size;
    }
    neighbours++; // N neighbours other than ourselves

    typedef KDTreeVectorOfVectorsAdaptor<std::vector<Vector3D>, double,
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
        neighbour_idx[idx].push_back(ret_indexes[i]);
      }
      idx++;
    }

    return neighbour_idx;
  }

  /// Remove outliers from Vector<Point3d> using Statistical Outlier algorithm
  ///
  /// @param points vector of points to filter
  /// @param neighbours Number of neighbours to consider. If less than 1 or
  /// greater than the number of points in the point cloud use all points
  /// @param outlier_margin Number of standard deviations
  /// @param verbose give verbose detail
  static void statistical_outlier_remover(std::vector<Vector3D> &points,
                                          size_t neighbours,
                                          double outlier_margin,
                                          bool verbose = false)
  {
    Timer("StatisticalOurtierRemover");
    std::vector<size_t> outliers =
        statistical_outlier_finder(points, neighbours, outlier_margin, verbose);
    std::vector<Vector3D> new_points;
    size_t k = 0;
    for (size_t i = 0; i < points.size(); i++)
    {
      if (k >= outliers.size() || i != outliers[k])
      {
        new_points.push_back(points[i]);
      }
      else
      {
        k++;
      }
    }
    points = new_points;
  }

  /// Remove outliers from PointCloud using Statistical Outlier algorithm
  ///
  /// @param point_cloud The point cloud to filter
  /// @param neighbours Number of neighbours to consider. If less than 1 or
  /// greater than the number of points in the point cloud use all points
  /// @param outlier_margin Number of standard deviations
  /// @param verbose give verbose detail
  static void statistical_outlier_remover(PointCloud &point_cloud,
                                          size_t neighbours,
                                          double outlier_margin,
                                          bool verbose = false)
  {
    // points to remove
    std::vector<size_t> outliers = statistical_outlier_finder(
        point_cloud.points, neighbours, outlier_margin, verbose);

    filter_point_cloud(point_cloud, outliers);
  }

  static void ransac_outlier_remover(PointCloud &point_cloud,
                                     double distance_threshold,
                                     size_t iterations = 100)
  {
    Timer("ransac_outlier_remover");
    std::vector<size_t> outliers = ransac_outlier_finder(
        point_cloud.points, distance_threshold, iterations);
    filter_point_cloud(point_cloud, outliers);
  }

  static void ransac_outlier_remover(std::vector<Vector3D> &points,
                                     double distance_threshold,
                                     size_t iterations = 100)
  {
    Timer("ransac_outlier_remover");

    auto outliers =
        ransac_outlier_finder(points, distance_threshold, iterations);
    if (outliers.size() == 0)
      return;
    std::vector<Vector3D> new_points;
    size_t k = 0;
    for (size_t i = 0; i < points.size(); i++)
    {
      if (k >= outliers.size() || i != outliers[k])
      {
        new_points.push_back(points[i]);
      }
      else
      {
        k++;
      }
    }
    points = new_points;
  }

  static std::vector<size_t>
  ransac_outlier_finder(std::vector<Vector3D> &points,
                        double distance_threshold,
                        size_t iterations = 100)
  {

    std::vector<size_t> outliers;
    if (points.size() < 9)
      return outliers;
    std::vector<size_t> best_outliers(points.size(), 0);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    std::default_random_engine generator(seed);
    std::uniform_int_distribution<size_t> distribution(0, points.size() - 1);
    auto rand_idx = std::bind(distribution, generator);
    Vector3D pt1, pt2, pt3;
    Vector3D v1, v2, v3;
    size_t idx1, idx2, idx3;
    double k;

    for (size_t i = 0; i < iterations; i++)
    {
      idx1 = rand_idx();
      idx2 = rand_idx();
      // idx1, 2 and 3 must be different
      while (idx2 == idx1)
      {
        idx2 = rand_idx();
      }
      idx3 = rand_idx();
      while (idx3 == idx1 || idx3 == idx2)
      {
        idx3 = rand_idx();
      }
      pt1 = points[idx1];
      pt2 = points[idx2];
      pt3 = points[idx3];
      v1 = Vector3D(pt2.x - pt1.x, pt2.y - pt1.y, pt2.z - pt1.z);
      v2 = Vector3D(pt3.x - pt1.x, pt3.y - pt1.y, pt3.z - pt1.z);
      v3 = v1.cross(v2);
      v3 /= v3.magnitude();
      if (isnan(v3.x)) // all three points are in a line
      {
        continue;
      }

      k = v3.dot(pt2);
      outliers.clear();
      for (size_t j = 0; j < points.size(); j++)
      {
        auto pt_plane_dist = std::abs(v3.dot(points[j]) - k) / v3.magnitude();
        if (pt_plane_dist > distance_threshold)
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
  static void filter_point_cloud(PointCloud &point_cloud,
                                 std::vector<size_t> pts_to_remove)
  {
    std::sort(pts_to_remove.begin(), pts_to_remove.end());
    std::vector<Vector3D> new_points;
    std::vector<Vector3D> new_normals;
    std::vector<Color> new_colors{};
    std::vector<uint8_t> new_classifications{};
    std::vector<uint16_t> new_intensities{};
    std::vector<uint8_t> new_scan_flags{};

    bool has_normals = point_cloud.normals.size() > 0;
    bool has_color = point_cloud.colors.size() > 0;
    bool has_class = point_cloud.classifications.size() > 0;
    bool has_intensity = point_cloud.intensities.size() > 0;
    bool has_scanflags = point_cloud.scan_flags.size() > 0;

    size_t k = 0;
    for (size_t i = 0; i < point_cloud.points.size(); i++)
    {
      if (k >= pts_to_remove.size() || i != pts_to_remove[k])
      {
        new_points.push_back(point_cloud.points[i]);

        if (has_normals)
          new_normals.push_back(point_cloud.normals[i]);
        if (has_color)
          new_colors.push_back(point_cloud.colors[i]);
        if (has_class)
          new_classifications.push_back(point_cloud.classifications[i]);
        if (has_intensity)
          new_intensities.push_back(point_cloud.intensities[i]);
        if (has_scanflags)
          new_scan_flags.push_back(point_cloud.scan_flags[i]);
      }
      else
      {
        k++;
      }
    }
    // Assign new to old
    point_cloud.points = new_points;

    if (has_normals)
      point_cloud.normals = new_normals;
    if (has_color)
      point_cloud.colors = new_colors;
    if (has_class)
      point_cloud.classifications = new_classifications;
    if (has_intensity)
      point_cloud.intensities = new_intensities;
    if (has_scanflags)
      point_cloud.scan_flags = new_scan_flags;
  }

  static std::pair<uint8_t, uint8_t> parse_scan_flag(uint8_t flag)
  {
    uint8_t return_number = flag & 7;
    uint8_t num_returns = (flag >> 3) & 7;
    return std::pair<uint8_t, uint8_t>(return_number, num_returns);
  }

  static uint8_t pack_scan_flag(uint8_t return_number, uint8_t num_returns)
  {
    return (return_number & 7) | ((num_returns & 7) << 3);
  }

  static PointCloud remove_vegetation(const PointCloud &point_cloud)
  {
    // Create copy of point cloud
    PointCloud _point_cloud{point_cloud};

    // Remove some points that might be vegetation
    std::vector<size_t> points_to_remove;
    bool has_scan_flags = true;
    if (_point_cloud.scan_flags.size() != _point_cloud.points.size())
    {
      warning("Scan flags not set, only using classification");
      has_scan_flags = false;
    }
    bool has_classification = false;
    if (_point_cloud.classifications.size() == _point_cloud.points.size())
      has_classification = true;

    if (!has_scan_flags && !has_classification)
    {
      warning("No scan flags or classification. No vegetation filtering");
      return _point_cloud;
    }

    for (size_t i = 0; i < _point_cloud.points.size(); i++)
    {
      auto scan_flag = std::pair<uint8_t, uint8_t>(0, 0);
      if (has_scan_flags)
        scan_flag = parse_scan_flag(_point_cloud.scan_flags[i]);

      uint8_t classification = 1;
      if (has_classification)
        classification = _point_cloud.classifications[i];
      // not last point of several

      if ((classification >= 3 && classification <= 5) || // classified as veg
          (classification < 2 &&                          // not classified
           (scan_flag.first != scan_flag.second))         // not last
      )
      {
        points_to_remove.push_back(i);
      }
    }
    filter_point_cloud(_point_cloud, points_to_remove);

    return _point_cloud;
  }

  static void estimate_normals_knn(PointCloud &point_cloud, size_t neighbours)
  {
    auto normals = estimate_normals_knn(point_cloud.points, neighbours);
    point_cloud.normals = normals;
  };

  static std::vector<Vector3D>
  estimate_normals_knn(std::vector<Vector3D> points, size_t neighbours)
  {
    std::vector<Vector3D> normals;
    auto neigbours_idx = knn_nearest_neighbours_idx(points, neighbours);
    size_t idx = 0;
    Eigen::RowVector3d dir(0.0, 0.0, 1.0);
    for (auto const &query_pt : points)
    {
      size_t found = neigbours_idx[idx].size();
      if (found < 3) // not enough neighbours to estimate normal
      {
        normals.push_back(Vector3D(0, 0, 0));
        idx++;
        continue;
      }

      Eigen::MatrixXd neighbors(found, 3);
      for (size_t i = 0; i < found; i++)
      {
        auto pt = points[neigbours_idx[idx][i]];
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
