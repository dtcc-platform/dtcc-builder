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

namespace VirtualCity
{

class HeightMapGenerator
{
public:

    // Generate height map using default method (nearest neighbor sampling)
    static void GenerateHeightMap(HeightMap& heightMap,
                                  const PointCloud& pointCloud)
    {
        GenerateHeightMapNearestNeighbor(heightMap, pointCloud);
    }

    // Generate height map using nearest neighbor sampling
    static void GenerateHeightMapNearestNeighbor(HeightMap& heightMap,
            const PointCloud& pointCloud)
    {
        std::cout << "HeightMapGenerator: Generating heightmap (nearest neighbor)..." << std::endl;

        // Initialize markers for closest points
        size_t numGridPoints = heightMap.GridData.size();
        std::vector<Point2D> closestPoints(numGridPoints);
        std::vector<bool> pointAdded(numGridPoints);
        std::fill(pointAdded.begin(), pointAdded.end(), false);

        // FIXME: Testing
        auto pp = *pointCloud.Points.begin();
        double xmin = pp.x;
        double xmax = pp.x;
        double ymin = pp.y;
        double ymax = pp.y;

        // Iterate over point cloud
        for (auto const & q3D : pointCloud.Points)
        {
            // Compute closest point in grid
            const Point2D q2D(q3D.x, q3D.y);
            const size_t i = heightMap.Coordinate2Index(q2D);

            // FIXME: Testing
            xmin = std::min(xmin, q2D.x);
            xmax = std::max(xmax, q2D.x);
            ymin = std::min(ymin, q2D.y);
            ymax = std::max(ymax, q2D.y);
            // std::cout << "q = " << q2D << std::endl;
            // std::cout << "i = " << i << std::endl;
            // std::cout << "p = " << heightMap.Index2Coordinate(i) << std::endl;
            // std::cout << std::endl;

            // Add if first point or closest
            if (!pointAdded[i])
            {
                heightMap.GridData[i] = q3D.z;
                pointAdded[i] = true;
            }
            else
            {
                const Point2D p2D = heightMap.Index2Coordinate(i);
                const Point2D& r2D = closestPoints[i];
                const double d2q = Geometry::SquaredDistance2D(p2D, q2D);
                const double d2r = Geometry::SquaredDistance2D(p2D, r2D);
                if (d2q < d2r)
                {
                    heightMap.GridData[i] = q3D.z;
                    closestPoints[i] = q2D;
                }
            }
        }


        std::cout << xmin << " " << ymin << " " << xmax << " " << ymax << std::endl;



        // Compute some stats
        size_t numAssigned = 0;
        double largestDistance = 0.0;
        for (size_t i = 0; i < numGridPoints; i++)
        {
            if (pointAdded[i])
            {
                numAssigned += 1;
                const Point2D& p2D = heightMap.Index2Coordinate(i);
                const Point2D& q2D = closestPoints[i];
                const double d2 = Geometry::SquaredDistance2D(p2D, q2D);
                if (d2 > largestDistance)
                    largestDistance = d2;
            }
        }
        largestDistance = std::sqrt(largestDistance);
        const double percentAssigned
            = 100.0 * static_cast<double>(numAssigned) / numGridPoints;
        std::cout << "HeightMapGenerator: "
                  << numGridPoints << " grid points" << std::endl;
        std::cout << "HeightMapGenerator: "
                  << numAssigned << " assigned points ("
                  << std::setprecision(3) <<percentAssigned
                  << "%)" << std::endl;
        std::cout << "HeightMapGenerator: largest distance = "
                  << largestDistance << std::endl;
    }

};

}

#endif
