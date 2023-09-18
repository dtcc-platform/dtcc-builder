// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_ELEVATION_BUILDER_H
#define DTCC_ELEVATION_BUILDER_H

#include <iomanip>
#include <iostream>
#include <stack>
#include <vector>

#include "Constants.h"
#include "Logging.h"
#include "Timer.h"
#include "VertexSmoother.h"
#include "model/GridField.h"
#include "model/PointCloud.h"
#include "model/Vector.h"

namespace DTCC_BUILDER
{

class ElevationBuilder
{
public:
  /// build digital elevation model (DEM) from point cloud.
  /// Only points matching the given classification(s) are used.
  /// If the classifications are empty, then all points are used.
  ///
  /// @param point_cloud Point cloud (unfiltered)
  /// @param classifications classifications to be considered
  /// @param resolution Resolution (grid size) of digital elevation model
  static GridField build_elevation(const PointCloud &point_cloud,
                                   const std::vector<int> &classifications,
                                   double resolution)
  {
    info("Buildinga digital elevation model from point cloud...");
    Timer timer("build_elevation");

    // Check that point cloud is not empty
    if (point_cloud.points.empty())
      error("empty point cloud");

    // Check that point cloud has classifications
    bool has_classification =
        (point_cloud.points.size() == point_cloud.classifications.size());
    if (!has_classification)
      warning("Missing classifications for point cloud, using all points");

    // print classifications
    if (has_classification && classifications.size() > 0)
    {
      std::string msg{"ElevationBuilder: Using classifications "};
      for (size_t i = 0; i < classifications.size(); i++)
      {
        msg += str(classifications[i]);
        if (i + 1 < classifications.size())
          msg += ", ";
      }
      info(msg);
    }
    else
    {
      info("Using all classifications");
    }

    // FIXME: Add function Grid::Init(bbox) that takes care of this
    // initialization.
    // FIXME: Use here and also in RandomizeElevation() below.

    // Initialize grid bounding box
    GridField dem;
    dem.grid.bounding_box = point_cloud.bounding_box;

    // Initialize grid data
    dem.grid.xsize =
        (dem.grid.bounding_box.Q.x - dem.grid.bounding_box.P.x) / resolution +
        1;
    dem.grid.ysize =
        (dem.grid.bounding_box.Q.y - dem.grid.bounding_box.P.y) / resolution +
        1;
    dem.values.resize(dem.grid.xsize * dem.grid.ysize);
    std::fill(dem.values.begin(), dem.values.end(), 0.0);
    dem.grid.xstep = (dem.grid.bounding_box.Q.x - dem.grid.bounding_box.P.x) /
                     (dem.grid.xsize - 1);
    dem.grid.ystep = (dem.grid.bounding_box.Q.y - dem.grid.bounding_box.P.y) /
                     (dem.grid.ysize - 1);

    // Compute mean raw elevation (used for skipping outliers)
    double mean_elevation_raw = 0.0;
    size_t num_inside = 0;
    for (auto const &p_3d : point_cloud.points)
    {
      // Get 2D point
      const Vector2D p_2d{p_3d.x, p_3d.y};

      // Skip if outside of domain
      if (!Geometry::bounding_box_contains_2d(dem.grid.bounding_box, p_2d))
        continue;

      // Sum up elevation
      mean_elevation_raw += p_3d.z;
      num_inside++;
    }
    mean_elevation_raw /= static_cast<double>(num_inside);

    // Initialize counters for number of points for local mean
    size_t num_grid_points = dem.values.size();
    std::vector<size_t> num_local_points(num_grid_points);
    std::fill(num_local_points.begin(), num_local_points.end(), 0);

    // Iterate over point cloud and sum up heights
    size_t num_outliers = 0;
    double mean_elevation = 0.0;
    num_inside = 0;
    std::vector<size_t> neighbor_indices;
    neighbor_indices.reserve(5);
    for (size_t i = 0; i < point_cloud.points.size(); i++)
    {
      // Get point and classification
      const Vector3D &p_3d{point_cloud.points[i]};
      uint8_t clf = 0;
      if (has_classification)
        clf = {point_cloud.classifications[i]};

      // Get 2D Point
      const Vector2D p_2d{p_3d.x, p_3d.y};

      // Check classification (accept all classifications if empty list)
      bool match = true;
      if (has_classification && classifications.size() > 0)
      {
        match = false;
        for (const auto c : classifications)
        {
          if (clf == c)
          {
            match = true;
            break;
          }
        }
      }

      // Skip if not matching classification
      if (!match)
        continue;

      // Skip if outside of domain
      if (!Geometry::bounding_box_contains_2d(dem.grid.bounding_box, p_2d))
        continue;

      // Skip if outlier
      if (p_3d.z - mean_elevation_raw >
          Constants::point_cloud_outlier_threshold)
      {
        num_outliers += 1;
        continue;
      }

      // Sum up elevation
      mean_elevation += p_3d.z;
      num_inside++;

      // Iterate over closest stencil (including center of stencil)
      neighbor_indices.clear();
      const size_t j = dem.grid.point_to_index(p_2d);
      neighbor_indices.push_back(j);
      dem.grid.index_to_boundary(j, neighbor_indices);
      for (size_t k : neighbor_indices)
      {
        dem.values[k] += p_3d.z;
        num_local_points[k] += 1;
      }
    }

    // Compute mean elevation
    mean_elevation /= static_cast<double>(num_inside);

    // Compute mean of elevations for each grid point
    std::vector<size_t> missing_indices;
    for (size_t i = 0; i < num_grid_points; i++)
    {
      if (num_local_points[i] > 0)
        dem.values[i] /= num_local_points[i];
      else
        missing_indices.push_back(i);
    }

    // Check that we have at least one point (very loose check)
    const size_t num_missing = missing_indices.size();
    if (num_missing == num_grid_points)
      throw std::runtime_error("No points inside height map domain.");

    // Reuse vector num_local_points to indicate which points have been
    // visited: 0 = empty, 1 = boundary, 2 = filled
    for (size_t i = 0; i < num_grid_points; i++)
      num_local_points[i] = (num_local_points[i] == 0 ? 0 : 2);

    // Create stack of boundary points neighboring unfilled regions by
    // examining the neighbors of all missing points. Note that we use
    // num_local_points to keep track of which boundary that have already
    // been added to the stack; only add neighbors that already contain
    // a value and only add neighbors that have not been added before.
    std::stack<size_t> boundary_indices;
    for (size_t i : missing_indices)
    {
      neighbor_indices.clear();
      dem.grid.index_to_boundary(i, neighbor_indices);
      for (size_t j : neighbor_indices)
      {
        if (num_local_points[j] == 2)
        {
          boundary_indices.push(j);
          num_local_points[j] = 1;
        }
      }
    }

    // Flood fill values until stack is empty
    size_t num_found = 0;
    while (!boundary_indices.empty())
    {
      // Get boundary index from top of stack
      const size_t i = boundary_indices.top();
      boundary_indices.pop();

      // Propagate values to neighbors and add neighbor to stack
      neighbor_indices.clear();
      dem.grid.index_to_boundary(i, neighbor_indices);
      for (size_t j : neighbor_indices)
      {
        if (num_local_points[j] == 0)
        {
          dem.values[j] = dem.values[i];
          boundary_indices.push(j);
          num_local_points[j] = 1;
          num_found++;
        }
      }
    }

    // Check that we found data for all grid points
    if (num_found != num_missing)
      throw std::runtime_error("Unable to find data for all grid points.");

    // print some stats
    const double percent_missing =
        100.0 * static_cast<double>(num_missing) / num_grid_points;
    info(str(num_outliers) + " outlier(s) ignored");
    info("mean elevation is " + str(mean_elevation, 4) + "m");
    info(str(num_grid_points) + " grid points");
    info(str(num_missing) + " missing grid points (" + str(percent_missing, 3) +
         "%)");
    //    std::cout << "ElevationBuilder: "
    //        << "Maximum search distance is " << maxStep << std::endl;

    // Test data for verifying orientation, bump in lower left corner
    // for (size_t i = 0; i < dem.values.size(); i++)
    // {
    //     Vector2D p = dem.Index2Coordinate(i);
    //     const double dx = dem.XMax - dem.XMin;
    //     const double dy = dem.YMax - dem.YMin;
    //     const double x = (p.x - dem.XMin) / dx;
    //     const double y = (p.y - dem.YMin) / dy;
    //     dem.values[i] = x * (1 - x) * (1 - x) * y * (1 - y) * (1 -
    //     y);
    // }

    return dem;
  }

  /// build a random elevation model. Used for benchmarking.
  ///
  /// @param dem The digital elevation model (DEM)
  /// @param bbox Bounding box of domain
  /// @param resolution Resolution (grid size) of digital elevation model
  static void randomize_elevationl(GridField &dem,
                                   const BoundingBox2D &bbox,
                                   double resolution)
  {
    info("Randomizing elevation model...");

    // Some hard-coded building dimensions
    const double H = 50.0;      // Maximum hill height
    const double W0 = 50.0;     // Minimum hill width
    const double W1 = 250.0;    // Maximum hill width
    const size_t num_hills = 30; // Number of hills (Gaussian bumps)

    // Initialize grid bounding box
    dem.grid.bounding_box = bbox;

    // Initialize grid data
    dem.grid.xsize =
        (dem.grid.bounding_box.Q.x - dem.grid.bounding_box.P.x) / resolution +
        1;
    dem.grid.ysize =
        (dem.grid.bounding_box.Q.y - dem.grid.bounding_box.P.y) / resolution +
        1;
    dem.values.resize(dem.grid.xsize * dem.grid.ysize);
    std::fill(dem.values.begin(), dem.values.end(), 0.0);
    dem.grid.xstep = (dem.grid.bounding_box.Q.x - dem.grid.bounding_box.P.x) /
                     (dem.grid.xsize - 1);
    dem.grid.ystep = (dem.grid.bounding_box.Q.y - dem.grid.bounding_box.P.y) /
                     (dem.grid.ysize - 1);

    // Randomize hills
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> h;
    std::vector<double> w;
    const double X0 = dem.grid.bounding_box.P.x;
    const double Y0 = dem.grid.bounding_box.P.y;
    const double dx = dem.grid.bounding_box.Q.x - dem.grid.bounding_box.P.x;
    const double dy = dem.grid.bounding_box.Q.y - dem.grid.bounding_box.P.y;
    for (size_t i = 0; i < num_hills; i++)
    {
      x.push_back(X0 + dx * Utils::random());
      y.push_back(Y0 + dy * Utils::random());
      h.push_back(H * Utils::random());
      w.push_back(W0 + (W1 - W0) * Utils::random());
    }

    // Set heights
    std::fill(dem.values.begin(), dem.values.end(), 0.0);
    for (size_t i = 0; i < dem.values.size(); i++)
    {
      const Vector2D p = dem.grid.index_to_point(i);
      for (size_t k = 0; k < num_hills; k++)
      {
        const double dx = p.x - x[k];
        const double dy = p.y - y[k];
        dem.values[i] += h[k] * exp(-0.5 * (dx * dx + dy * dy) / (w[k] * w[k]));
      }
    }
  }
};

} // namespace DTCC_BUILDER

#endif
