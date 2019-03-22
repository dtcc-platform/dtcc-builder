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

        // Iterate over footprints
        for (auto const & polygon : polygons)
        {
            // Compute center and radius of footprint
            const Point2D center = polygon.Center();
            const double radius = polygon.Radius(center);

            // Add points for sampling height
            std::vector<Point2D> samplePoints;
            const double a = 0.5*radius;
            samplePoints.push_back(center);
            samplePoints.push_back(Point2D(center.x + a, center.y));
            samplePoints.push_back(Point2D(center.x - a, center.y));
            samplePoints.push_back(Point2D(center.x, center.y + a));
            samplePoints.push_back(Point2D(center.x, center.y - a));

            // Compute mean height at points inside footprints
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
                std::cout << "CityModelGenerator: Skipping building, no sample points inside building footprint" << std::endl;
                continue;
            }

            // Create building
            Building building;
            for (auto const & p : polygon.Points)
                building.Footprint.Points.push_back(p);
            building.Height = z / numInside;

            std::cout << building.Height << std::endl;

            // Add building
            cityModel.Buildings.push_back(building);
        }
    }

};

}

#endif
