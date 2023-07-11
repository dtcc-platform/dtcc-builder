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
  /// Build digital elevation model (DEM) from point cloud.
  /// Only points matching the given classification(s) are used.
  /// If the classifications are empty, then all points are used.
  ///
  /// @param pointCloud Point cloud (unfiltered)
  /// @param classifications Classifications to be considered
  /// @param resolution Resolution (grid size) of digital elevation model
  static GridField build_elevation(const PointCloud &pointCloud,
                                   const std::vector<int> &classifications,
                                   double resolution)
  {
    info("Buildinga digital elevation model from point cloud...");
    Timer timer("build_elevation");

    // Check that point cloud is not empty
    if (pointCloud.Points.empty())
      error("Empty point cloud");

    // Check that point cloud has classifications
    bool has_classification =
        (pointCloud.Points.size() == pointCloud.Classifications.size());
    if (!has_classification)
      warning("Missing classifications for point cloud, using all points");

    // Print classifications
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
    dem.grid.BoundingBox = pointCloud.BoundingBox;

    // Initialize grid data
    dem.grid.XSize =
        (dem.grid.BoundingBox.Q.x - dem.grid.BoundingBox.P.x) / resolution + 1;
    dem.grid.YSize =
        (dem.grid.BoundingBox.Q.y - dem.grid.BoundingBox.P.y) / resolution + 1;
    dem.Values.resize(dem.grid.XSize * dem.grid.YSize);
    std::fill(dem.Values.begin(), dem.Values.end(), 0.0);
    dem.grid.XStep = (dem.grid.BoundingBox.Q.x - dem.grid.BoundingBox.P.x) /
                     (dem.grid.XSize - 1);
    dem.grid.YStep = (dem.grid.BoundingBox.Q.y - dem.grid.BoundingBox.P.y) /
                     (dem.grid.YSize - 1);

    // Compute mean raw elevation (used for skipping outliers)
    double meanElevationRaw = 0.0;
    size_t numInside = 0;
    for (auto const &p3D : pointCloud.Points)
    {
      // Get 2D point
      const Point2D p2D{p3D.x, p3D.y};

      // Skip if outside of domain
      if (!Geometry::BoundingBoxContains2D(dem.grid.BoundingBox, p2D))
        continue;

      // Sum up elevation
      meanElevationRaw += p3D.z;
      numInside++;
    }
    meanElevationRaw /= static_cast<double>(numInside);

    // Initialize counters for number of points for local mean
    size_t numGridPoints = dem.Values.size();
    std::vector<size_t> numLocalPoints(numGridPoints);
    std::fill(numLocalPoints.begin(), numLocalPoints.end(), 0);

    // Iterate over point cloud and sum up heights
    size_t numOutliers = 0;
    double meanElevation = 0.0;
    numInside = 0;
    std::vector<size_t> neighborIndices;
    neighborIndices.reserve(5);
    for (size_t i = 0; i < pointCloud.Points.size(); i++)
    {
      // Get point and classification
      const Point3D &p3D{pointCloud.Points[i]};
      uint8_t clf = 0;
      if (has_classification)
        clf = {pointCloud.Classifications[i]};

      // Get 2D Point
      const Point2D p2D{p3D.x, p3D.y};

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
      if (!Geometry::BoundingBoxContains2D(dem.grid.BoundingBox, p2D))
        continue;

      // Skip if outlier
      if (p3D.z - meanElevationRaw > Constants::PointCloudOutlierThreshold)
      {
        numOutliers += 1;
        continue;
      }

      // Sum up elevation
      meanElevation += p3D.z;
      numInside++;

      // Iterate over closest stencil (including center of stencil)
      neighborIndices.clear();
      const size_t j = dem.grid.Point2Index(p2D);
      neighborIndices.push_back(j);
      dem.grid.Index2Boundary(j, neighborIndices);
      for (size_t k : neighborIndices)
      {
        dem.Values[k] += p3D.z;
        numLocalPoints[k] += 1;
      }
    }

    // Compute mean elevation
    meanElevation /= static_cast<double>(numInside);

    // Compute mean of elevations for each grid point
    std::vector<size_t> missingIndices;
    for (size_t i = 0; i < numGridPoints; i++)
    {
      if (numLocalPoints[i] > 0)
        dem.Values[i] /= numLocalPoints[i];
      else
        missingIndices.push_back(i);
    }

    // Check that we have at least one point (very loose check)
    const size_t numMissing = missingIndices.size();
    if (numMissing == numGridPoints)
      throw std::runtime_error("No points inside height map domain.");

    // Reuse vector numLocalPoints to indicate which points have been
    // visited: 0 = empty, 1 = boundary, 2 = filled
    for (size_t i = 0; i < numGridPoints; i++)
      numLocalPoints[i] = (numLocalPoints[i] == 0 ? 0 : 2);

    // Create stack of boundary points neighboring unfilled regions by
    // examining the neighbors of all missing points. Note that we use
    // numLocalPoints to keep track of which boundary that have already
    // been added to the stack; only add neighbors that already contain
    // a value and only add neighbors that have not been added before.
    std::stack<size_t> boundaryIndices;
    for (size_t i : missingIndices)
    {
      neighborIndices.clear();
      dem.grid.Index2Boundary(i, neighborIndices);
      for (size_t j : neighborIndices)
      {
        if (numLocalPoints[j] == 2)
        {
          boundaryIndices.push(j);
          numLocalPoints[j] = 1;
        }
      }
    }

    // Flood fill values until stack is empty
    size_t numFound = 0;
    while (!boundaryIndices.empty())
    {
      // Get boundary index from top of stack
      const size_t i = boundaryIndices.top();
      boundaryIndices.pop();

      // Propagate values to neighbors and add neighbor to stack
      neighborIndices.clear();
      dem.grid.Index2Boundary(i, neighborIndices);
      for (size_t j : neighborIndices)
      {
        if (numLocalPoints[j] == 0)
        {
          dem.Values[j] = dem.Values[i];
          boundaryIndices.push(j);
          numLocalPoints[j] = 1;
          numFound++;
        }
      }
    }

    // Check that we found data for all grid points
    if (numFound != numMissing)
      throw std::runtime_error("Unable to find data for all grid points.");

    // Print some stats
    const double percentMissing =
        100.0 * static_cast<double>(numMissing) / numGridPoints;
    info(str(numOutliers) + " outlier(s) ignored");
    info("Mean elevation is " + str(meanElevation, 4) + "m");
    info(str(numGridPoints) + " grid points");
    info(str(numMissing) + " missing grid points (" + str(percentMissing, 3) +
         "%)");
    //    std::cout << "ElevationBuilder: "
    //        << "Maximum search distance is " << maxStep << std::endl;

    // Test data for verifying orientation, bump in lower left corner
    // for (size_t i = 0; i < dem.Values.size(); i++)
    // {
    //     Vector2D p = dem.Index2Coordinate(i);
    //     const double dx = dem.XMax - dem.XMin;
    //     const double dy = dem.YMax - dem.YMin;
    //     const double x = (p.x - dem.XMin) / dx;
    //     const double y = (p.y - dem.YMin) / dy;
    //     dem.Values[i] = x * (1 - x) * (1 - x) * y * (1 - y) * (1 -
    //     y);
    // }

    return dem;
  }

  /// Build a random elevation model. Used for benchmarking.
  ///
  /// @param dem The digital elevation model (DEM)
  /// @param bbox Bounding box of domain
  /// @param resolution Resolution (grid size) of digital elevation model
  static void RandomizeElevationl(GridField &dem,
                                  const BoundingBox2D &bbox,
                                  double resolution)
  {
    info("Randomizing elevation model...");

    // Some hard-coded building dimensions
    const double H = 50.0;      // Maximum hill height
    const double W0 = 50.0;     // Minimum hill width
    const double W1 = 250.0;    // Maximum hill width
    const size_t numHills = 30; // Number of hills (Gaussian bumps)

    // Initialize grid bounding box
    dem.grid.BoundingBox = bbox;

    // Initialize grid data
    dem.grid.XSize =
        (dem.grid.BoundingBox.Q.x - dem.grid.BoundingBox.P.x) / resolution + 1;
    dem.grid.YSize =
        (dem.grid.BoundingBox.Q.y - dem.grid.BoundingBox.P.y) / resolution + 1;
    dem.Values.resize(dem.grid.XSize * dem.grid.YSize);
    std::fill(dem.Values.begin(), dem.Values.end(), 0.0);
    dem.grid.XStep = (dem.grid.BoundingBox.Q.x - dem.grid.BoundingBox.P.x) /
                     (dem.grid.XSize - 1);
    dem.grid.YStep = (dem.grid.BoundingBox.Q.y - dem.grid.BoundingBox.P.y) /
                     (dem.grid.YSize - 1);

    // Randomize hills
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> h;
    std::vector<double> w;
    const double X0 = dem.grid.BoundingBox.P.x;
    const double Y0 = dem.grid.BoundingBox.P.y;
    const double dX = dem.grid.BoundingBox.Q.x - dem.grid.BoundingBox.P.x;
    const double dY = dem.grid.BoundingBox.Q.y - dem.grid.BoundingBox.P.y;
    for (size_t i = 0; i < numHills; i++)
    {
      x.push_back(X0 + dX * Utils::Random());
      y.push_back(Y0 + dY * Utils::Random());
      h.push_back(H * Utils::Random());
      w.push_back(W0 + (W1 - W0) * Utils::Random());
    }

    // Set heights
    std::fill(dem.Values.begin(), dem.Values.end(), 0.0);
    for (size_t i = 0; i < dem.Values.size(); i++)
    {
      const Point2D p = dem.grid.Index2Point(i);
      for (size_t k = 0; k < numHills; k++)
      {
        const double dx = p.x - x[k];
        const double dy = p.y - y[k];
        dem.Values[i] += h[k] * exp(-0.5 * (dx * dx + dy * dy) / (w[k] * w[k]));
      }
    }
  }
};

} // namespace DTCC_BUILDER

#endif
