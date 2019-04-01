// City model generation from building footprints and height map.
// Copyright (C) 2019 Anders Logg.

#ifndef VC_CITY_MODEL_GENERATOR_H
#define VC_CITY_MODEL_GENERATOR_H

#include <iostream>
#include <vector>

#include "Point.h"
#include "Polygon.h"
#include "CityModel.h"

namespace VirtualCity
{

class CityModelGenerator
{
public:

    // Generate city model from building footprints and height map
    static void GenerateCityModel(CityModel& cityModel,
                                  const std::vector<Polygon>& polygons,
                                  const HeightMap& heightMap)
    {
        std::cout << "CityModelGenerator: Generating city model..."
                  << std::endl;

        // Array of unique points
        std::vector<Point2D> uniquePoints;

        // Iterate over footprints
        for (auto const & polygon : polygons)
        {
            // Check that the polygon is closed
            const size_t numPoints = polygon.Points.size();
            const double d = Geometry::Distance2D(polygon.Points[0],
                                                  polygon.Points[numPoints - 1]);
            if (d > Parameters::Epsilon)
            {
                std::cout << "CityModelGenerator: Skipping building, expecting polygon to be closed." << std::endl;
                continue;
            }

            // Check for intersecting polygons
            if (Intersects(polygon.Points, uniquePoints))
            {
                std::cout << "CityModelGenerator: Skipping building, too close to existing building." << std::endl;
                continue;
            }
            else
            {
                for (auto const & p : polygon.Points)
                    uniquePoints.push_back(p);
            }

            // Compute center and radius of footprint
            const Point2D center = polygon.Center();
            const double radius = polygon.Radius(center);

            // Add points for sampling height
            std::vector<Point2D> samplePoints;
            const double a = 0.5 * radius;
            samplePoints.push_back(center);
            samplePoints.push_back(Point2D(center.x + a, center.y));
            samplePoints.push_back(Point2D(center.x - a, center.y));
            samplePoints.push_back(Point2D(center.x, center.y + a));
            samplePoints.push_back(Point2D(center.x, center.y - a));

            // Compute mean height at points inside footprint
            double z = 0.0;
            size_t numInside = 0;
            for (auto const & p : samplePoints)
            {
                if (polygon.Contains(p))
                {
                    z += heightMap(p);
                    numInside += 1;
                }
            }

            // Check if we got at least one point
            if (numInside == 0)
            {
                std::cout << "CityModelGenerator: Skipping building, no sample points inside building footprint." << std::endl;
                continue;
            }

            // Set building height
            Building building;
            building.Height = z / numInside;
            //std::cout << "Height = " << building.Height << std::endl;

            // Set building footprint (skip last point and any duplicates)
            for (size_t i = 0; i < numPoints - 1; i++)
            {
                // Check distance to previous point
                if (i > 0 && Geometry::Distance2D(polygon.Points[i],
                                                  polygon.Points[i - 1])
                        < Parameters::FootprintDuplicateThreshold)
                {
                    std::cout << "CityModelGenerator: Skipping duplicate point in footprint." << std::endl;
                    continue;
                }

                // Add point
                building.Footprint.Points.push_back(polygon.Points[i]);
            }

            // Add building
            cityModel.Buildings.push_back(building);
        }
    }

private:

    // Check if polygon intersects (or is very close to an existing point)
    static bool Intersects(const std::vector<Point2D>& polygon,
                           const std::vector<Point2D>& uniquePoints)
    {
        double eps2 = Parameters::FootprintDuplicateThreshold;
        eps2 = eps2 * eps2;
        for (auto const & p : polygon)
        {
            for (auto const & q : uniquePoints)
            {
                if (Geometry::SquaredDistance2D(p, q) < eps2)
                    return true;
            }
        }
        return false;
    }

};

}

#endif
