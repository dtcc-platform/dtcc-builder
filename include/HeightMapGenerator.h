// Height map generation from point cloud (LiDAR) data.
// Copyright (C) 2019 Anders Logg.

#ifndef DTCC_HEIGHT_MAP_GENERATOR_H
#define DTCC_HEIGHT_MAP_GENERATOR_H

#include <iomanip>
#include <iostream>
#include <vector>

#include "Geometry.h"
#include "HeightMap.h"
#include "Parameters.h"
#include "Point.h"
#include "PointCloud.h"

namespace DTCC
{

class HeightMapGenerator
{
public:
  // Generate height map from point cloud
  static void GenerateHeightMap(HeightMap &heightMap,
                                const PointCloud &pointCloud,
                                double x0,
                                double y0,
                                double xMin,
                                double yMin,
                                double xMax,
                                double yMax,
                                double heightMapResolution)
  {
    std::cout << "HeightMapGenerator: Generating heightmap from point cloud..."
              << std::endl;

    // Shortcut
    HeightMap &hm = heightMap;

    // Initialize grid dimensions
    hm.XMin = xMin;
    hm.YMin = yMin;
    hm.XMax = xMax;
    hm.YMax = yMax;

    // Initialize grid data
    hm.XSize = (hm.XMax - hm.XMin) / heightMapResolution + 1;
    hm.YSize = (hm.YMax - hm.YMin) / heightMapResolution + 1;
    hm.GridData.resize(hm.XSize * hm.YSize);
    std::fill(hm.GridData.begin(), hm.GridData.end(), 0.0);
    hm.XStep = (hm.XMax - hm.XMin) / (hm.XSize - 1);
    hm.YStep = (hm.YMax - hm.YMin) / (hm.YSize - 1);

    std::cout << "HeightMapGenerator: Computing mean elevation" << std::endl;

    // Compute mean raw elevation (used for skipping outliers)
    double meanElevationRaw = 0.0;
    for (auto const &q3D : pointCloud.Points)
      meanElevationRaw += q3D.z;
    meanElevationRaw /= pointCloud.Points.size();

    // Initialize counters for number of points for local mean
    size_t numGridPoints = hm.GridData.size();
    std::vector<size_t> numLocalPoints(numGridPoints);
    std::fill(numLocalPoints.begin(), numLocalPoints.end(), 0);

    std::cout << "HeightMapGenerator: Extracting point cloud data" << std::endl;

    // Iterate over point cloud and sum up heights
    size_t numOutliers = 0;
    double meanElevation = 0.0;
    for (auto const &q3D : pointCloud.Points)
    {
      // Ignore outliers
      if (q3D.z - meanElevationRaw > Parameters::PointCloudOutlierThreshold)
      {
        numOutliers += 1;
        continue;
      }

      // Get 2D point and subtract origin
      const Point2D q2D(q3D.x - x0, q3D.y - y0);

      // Recompute mean elevation (excluding outliers)
      meanElevation += q3D.z;

      // Iterate over neighbors in grid
      for (size_t i : hm.Coordinate2Indices(q2D))
      {
        // Compute distance to grid point
        const Point2D p2D = hm.Index2Coordinate(i);
        const double dx = std::abs(p2D.x - q2D.x);
        const double dy = std::abs(p2D.y - q2D.y);

        // Note: The test below should (almost) always be true
        // with current threshold but we might want to use a
        // tighter threshold.

        // Add if closer than threshold
        if (dx < hm.XStep && dy < hm.YStep)
        {
          hm.GridData[i] += q3D.z;
          numLocalPoints[i] += 1;
        }
      }
    }

    // Compute mean elevation
    meanElevation /= pointCloud.Points.size() - numOutliers;

    std::cout << "HeightMapGenerator: Computing local mean elevation"
              << std::endl;

    // Compute mean of elevations for each grid point
    std::vector<size_t> missingIndices;
    for (size_t i = 0; i < numGridPoints; i++)
    {
      if (numLocalPoints[i] > 0)
        hm.GridData[i] /= numLocalPoints[i];
      else
        missingIndices.push_back(i);
    }

    // Check that we have at least one point (very loose check)
    const size_t numMissing = missingIndices.size();
    if (numMissing == numGridPoints)
      throw std::runtime_error("No points inside height map domain.");

    // Note: We fill in missing point by searching for the
    // closest existing value around each missing grid point.
    // It might be more efficient to do a flood fill.

    std::cout << "HeightMapGenerator: Filling in missing grid points ("
              << numMissing << "/" << numGridPoints << ")" << std::endl;

    // Fill in data for missing grid points
    const size_t maxDiameter = std::max(hm.XSize, hm.YSize);
    size_t numFound = 0;
    size_t maxStep = 0;
    for (size_t i : missingIndices)
    {
      // Iterate over larger and larger boundaries
      for (size_t step = 1; step < maxDiameter; step++)
      {
        // Compute maximum step
        if (step > maxStep)
          maxStep = step;

        // Search grid points at distance step
        bool found = false;
        for (size_t j : hm.Index2Boundary(i, step))
        {
          if (numLocalPoints[j] > 0)
          {
            hm.GridData[i] = hm.GridData[j];
            numFound += 1;
            found = true;
            break;
          }
        }

        // Point found
        if (found)
          break;
      }
    }

    // Check that we found data for all grid points
    if (numFound != numMissing)
      throw std::runtime_error("Unable to find data for all grid points.");

    // Print some stats
    const double percentMissing =
        100.0 * static_cast<double>(numMissing) / numGridPoints;
    std::cout << "HeightMapGenerator: " << numOutliers << " outliers ignored"
              << std::endl;
    std::cout << "HeightMapGenerator: Mean elevation is "
              << std::setprecision(4) << meanElevation << "m" << std::endl;
    std::cout << "HeightMapGenerator: " << numGridPoints << " grid points"
              << std::endl;
    std::cout << "HeightMapGenerator: " << numMissing
              << " missing grid points (" << std::setprecision(3)
              << percentMissing << "%)" << std::endl;
    std::cout << "HeightMapGenerator: "
              << "Maximum search distance is " << maxStep << std::endl;

    // Test data for verifying orientation, bump in lower left corner
    // for (size_t i = 0; i < heightMap.GridData.size(); i++)
    // {
    //     Point2D p = heightMap.Index2Coordinate(i);
    //     const double dx = heightMap.XMax - heightMap.XMin;
    //     const double dy = heightMap.YMax - heightMap.YMin;
    //     const double x = (p.x - heightMap.XMin) / dx;
    //     const double y = (p.y - heightMap.YMin) / dy;
    //     heightMap.GridData[i] = x * (1 - x) * (1 - x) * y * (1 - y) * (1 -
    //     y);
    // }
  }
};

} // namespace DTCC

#endif
