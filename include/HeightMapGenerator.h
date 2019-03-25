// Height map generation from point cloud (LiDAR) data.
// Copyright (C) 2019 Anders Logg.

#ifndef VC_HEIGHT_MAP_GENERATOR_H
#define VC_HEIGHT_MAP_GENERATOR_H

#include <iostream>
#include <iomanip>
#include <vector>


#include "Point.h"
#include "Geometry.h"
#include "HeightMap.h"
#include "PointCloud.h"
#include "Parameters.h"

namespace VirtualCity
{

class HeightMapGenerator
{
public:

    // Generate height map from point cloud
    static void GenerateHeightMap(HeightMap& heightMap,
                                  const PointCloud& pointCloud)
    {
        std::cout << "HeightMapGenerator: Generating heightmap from point cloud..." << std::endl;

        // Compute mean height
        std::cout << "HeightMapGenerator: Computing mean elevation" << std::endl;
        double meanElevationRaw = 0.0;
        for (auto const & q3D : pointCloud.Points)
            meanElevationRaw += q3D.z;
        meanElevationRaw /= pointCloud.Points.size();

        // Initialize counters for number of points for local mean
        std::vector<double>& gridData = heightMap.GridData;
        size_t numGridPoints = gridData.size();
        std::vector<size_t> numLocalPoints(numGridPoints);
        std::fill(numLocalPoints.begin(), numLocalPoints.end(), 0);

        // Reset grid data
        std::fill(gridData.begin(), gridData.end(), 0.0);

        std::cout << "HeightMapGenerator: Extracting point cloud data"
                  << std::endl;

        // Iterate over point cloud
        size_t numOutliers = 0;
        double meanElevation = 0.0;
        for (auto const & q3D : pointCloud.Points)
        {
            // Ignore outliers
            if (q3D.z - meanElevationRaw > Parameters::PointCloudOutlierThreshold)
            {
                numOutliers += 1;
                continue;
            }

            // Get 2D point
            const Point2D q2D(q3D.x, q3D.y);

            // Recompute mean elevation (excluding outliers)
            meanElevation += q3D.z;

            // Iterate over neighbors in grid
            for (size_t i : heightMap.Coordinate2Indices(q2D))
            {
                // Compute distance to grid point
                const Point2D p2D = heightMap.Index2Coordinate(i);
                const double dx = std::abs(p2D.x - q2D.x);
                const double dy = std::abs(p2D.y - q2D.y);

                // Note: The test below should (almost) always be true
                // with current threshold but we might want to use a
                // tighter threshold.

                // Add if closer than threshold
                if (dx < heightMap.XStep && dy < heightMap.YStep)
                {
                    gridData[i] += q3D.z;
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
                gridData[i] /= numLocalPoints[i];
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

        std::cout << "HeightMapGenerator: Filling in missing grid points"
                  << std::endl;

        // Fill in data for missing grid points
        const size_t maxDiameter = std::max(heightMap.XSize,
                                            heightMap.YSize);
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
                for (size_t j : heightMap.Index2Boundary(i, step))
                {
                    if (numLocalPoints[j] > 0)
                    {
                        gridData[i] = gridData[j];
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
        const double percentMissing
            = 100.0 * static_cast<double>(numMissing) / numGridPoints;
        std::cout << "HeightMapGenerator: "
                  << numOutliers << " outliers ignored" << std::endl;
        std::cout << "HeightMapGenerator: Mean elevation is "
                  << std::setprecision(4) << meanElevation
                  << "m" << std::endl;
        std::cout << "HeightMapGenerator: "
                  << numGridPoints << " grid points" << std::endl;
        std::cout << "HeightMapGenerator: "
                  << numMissing << " missing grid points ("
                  << std::setprecision(3) << percentMissing
                  << "%)" << std::endl;
        std::cout << "HeightMapGenerator: "
                  << "Maximum search distance is "
                  << maxStep << std::endl;

        // Test data for verifying orientation, bump in lower left corner
        // for (size_t i = 0; i < heightMap.GridData.size(); i++)
        // {
        //     Point2D p = heightMap.Index2Coordinate(i);
        //     const double dx = heightMap.XMax - heightMap.XMin;
        //     const double dy = heightMap.YMax - heightMap.YMin;
        //     const double x = (p.x - heightMap.XMin) / dx;
        //     const double y = (p.y - heightMap.YMin) / dy;
        //     heightMap.GridData[i] = x * (1 - x) * (1 - x) * y * (1 - y) * (1 - y);
        // }

    }

};

}

#endif
