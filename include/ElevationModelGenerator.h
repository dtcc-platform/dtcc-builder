// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_ELEVATION_MODEL_GENERATOR_H
#define DTCC_ELEVATION_MODEL_GENERATOR_H

#include <iomanip>
#include <iostream>
#include <stack>
#include <vector>

#include "GridField.h"
#include "Logging.h"
#include "Parameters.h"
#include "PointCloud.h"
#include "Timer.h"
#include "Vector.h"

namespace DTCC
{

class ElevationModelGenerator
{
public:
  // Generate digital elevation model (DEM) from point cloud
  static void GenerateElevationModel(GridField2D &dem,
                                     const PointCloud &pointCloud,
                                     double resolution,
                                     std::string type = "DEM")
  {
    Info("ElevationModelGenerator: Generating digital elevation model (" +
         type + ") from point cloud...");
    Timer("GenerateElevationModel");

    // Check for empty data
    if (pointCloud.Points.empty())
      Error("ElevationModelGenerator: Empty point cloud");

    // Initialize grid bounding box
    dem.Grid.BoundingBox = pointCloud.BoundingBox;

    // Initialize grid data
    dem.Grid.XSize =
        (dem.Grid.BoundingBox.Q.x - dem.Grid.BoundingBox.P.x) / resolution + 1;
    dem.Grid.YSize =
        (dem.Grid.BoundingBox.Q.y - dem.Grid.BoundingBox.P.y) / resolution + 1;
    dem.Values.resize(dem.Grid.XSize * dem.Grid.YSize);
    std::fill(dem.Values.begin(), dem.Values.end(), 0.0);
    dem.Grid.XStep = (dem.Grid.BoundingBox.Q.x - dem.Grid.BoundingBox.P.x) /
                     (dem.Grid.XSize - 1);
    dem.Grid.YStep = (dem.Grid.BoundingBox.Q.y - dem.Grid.BoundingBox.P.y) /
                     (dem.Grid.YSize - 1);

    Progress("ElevationModelGenerator: Computing mean elevation");

    // Compute mean raw elevation (used for skipping outliers)
    double meanElevationRaw = 0.0;
    size_t numInside = 0;
    for (auto const &q3D : pointCloud.Points)
    {
      // Get 2D point
      const Vector2D q2D(q3D.x, q3D.y);

      // Skip if outside of domain
      if (!Geometry::BoundingBoxContains2D(dem.Grid.BoundingBox, q2D))
        continue;

      // Sum up elevation
      meanElevationRaw += q3D.z;
      numInside++;
    }
    meanElevationRaw /= static_cast<double>(numInside);

    // Initialize counters for number of points for local mean
    size_t numGridPoints = dem.Values.size();
    std::vector<size_t> numLocalPoints(numGridPoints);
    std::fill(numLocalPoints.begin(), numLocalPoints.end(), 0);

    Progress("ElevationModelGenerator: Extracting point cloud data");

    // Iterate over point cloud and sum up heights
    size_t numOutliers = 0;
    double meanElevation = 0.0;
    numInside = 0;
    std::vector<size_t> neighborIndices;
    neighborIndices.reserve(5);
    for (auto const &q3D : pointCloud.Points)
    {
      // Get 2D point
      const Vector2D q2D(q3D.x, q3D.y);

      // Skip if outside of domain
      if (!Geometry::BoundingBoxContains2D(dem.Grid.BoundingBox, q2D))
        continue;

      // Ignore outliers
      if (q3D.z - meanElevationRaw > Parameters::PointCloudOutlierThreshold)
      {
        numOutliers += 1;
        continue;
      }

      // Sum up elevation
      meanElevation += q3D.z;
      numInside++;

      // Iterate over closest stencil (including center of stencil)
      neighborIndices.clear();
      const size_t i = dem.Grid.Point2Index(q2D);
      neighborIndices.push_back(i);
      dem.Grid.Index2Boundary(i, neighborIndices);
      for (size_t j : neighborIndices)
      {
        dem.Values[j] += q3D.z;
        numLocalPoints[j] += 1;
      }
    }

    // Compute mean elevation
    meanElevation /= static_cast<double>(numInside);

    Progress("ElevationModelGenerator: Computing local mean elevation");

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

    Progress("ElevationModelGenerator: Filling in missing grid points (" +
             str(numMissing) + "/" + str(numGridPoints) + ")");

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
      dem.Grid.Index2Boundary(i, neighborIndices);
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
      dem.Grid.Index2Boundary(i, neighborIndices);
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
    Info("ElevationModelGenerator: " + str(numOutliers) + " outliers ignored");
    Info("ElevationModelGenerator: Mean elevation is " + str(meanElevation, 4) +
         "m");
    Info("ElevationModelGenerator: " + str(numGridPoints) + " grid points");
    Info("ElevationModelGenerator: " + str(numMissing) +
         " missing grid points (" + str(percentMissing, 3) + "%)");
    //    std::cout << "ElevationModelGenerator: "
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
  }
};

} // namespace DTCC

#endif
