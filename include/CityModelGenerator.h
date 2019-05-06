// City model generation from building footprints and height map.
// Copyright (C) 2019 Anders Logg.

#ifndef VC_CITY_MODEL_GENERATOR_H
#define VC_CITY_MODEL_GENERATOR_H

#include <iostream>
#include <vector>
#include <queue>

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
                                  const HeightMap& heightMap,
                                  double minimalBuildingDistance)
    {
        std::cout << "CityModelGenerator: Generating city model..."
                  << std::endl;

        // Copy polygon data
        std::vector<Polygon> _polygons = polygons;

        // Compute closed polygons
        _polygons = ComputeClosedPolygons(_polygons);

        // Compute counter-clockwise oriented polygons
        _polygons = ComputeOrientedPolygons(_polygons);

        // Compute merged polygons
        _polygons = ComputeMergedPolygons(_polygons, minimalBuildingDistance);

        // Add buildings
        for (auto const & polygon : _polygons)
        {
            Building building;
            building.Footprint = polygon;
            cityModel.Buildings.push_back(building);
        }

    }


    /*
            // Array of unique points
            std::vector<Point2D> uniquePoints;

            // Iterate over polygons
            for (auto const & polygon : polygons)
            {


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


    */

private:

    // Compute closed polygons (removing any duplicate vertices)
    static std::vector<Polygon>
    ComputeClosedPolygons(const std::vector<Polygon>& polygons)
    {
        std::vector<Polygon> closedPolygons;
        size_t numSkipped = 0;

        // Iterate over polygons
        for (auto const & polygon : polygons)
        {
            // Compute distance between first and last point
            const size_t numPoints = polygon.Points.size();
            const double d = Geometry::Distance2D(polygon.Points[0],
                                                  polygon.Points[numPoints - 1]);

            // Skip if not closed
            if (d > Parameters::Epsilon)
            {
                numSkipped++;
                continue;
            }

            // Add polygon but include last (duplicate) point
            Polygon closedPolygon;
            for (size_t i = 0; i < polygon.Points.size() - 1; i++)
                closedPolygon.Points.push_back(polygon.Points[i]);
            closedPolygons.push_back(closedPolygon);
        }

        std::cout << "CityModelGenerator: Skipped " << numSkipped
                  << " building(s); expecting polygons to be closed."
                  << std::endl;

        return closedPolygons;
    }

    // Compute counter-clockwise oriented polygons. It is assumed that the
    // input polygons are closed without duplicate vertices.
    static std::vector<Polygon>
    ComputeOrientedPolygons(const std::vector<Polygon>& polygons)
    {
        std::vector<Polygon> orientedPolygons;
        size_t numReversed = 0;

        // Iterate over polygons
        for (auto const & polygon : polygons)
        {
            // Copy polygon
            Polygon orientedPolygon = polygon;

            // Reverse polygon if not counter-clockwise
            if (Geometry::PolygonOrientation2D(orientedPolygon) != 0)
            {
                numReversed++;
                std::reverse(orientedPolygon.Points.begin(),
                             orientedPolygon.Points.end());
            }

            // Add polygon
            orientedPolygons.push_back(orientedPolygon);
        }

        std::cout << "CityModelGenerator: Reversed " << numReversed
                  << " polygon(s) out of " << polygons.size() << "."
                  << std::endl;

        return orientedPolygons;
    }

    // Compute merged polygons. It is assumed that the input polygons
    // are closed without duplicate vertices and oriented.
    static std::vector<Polygon>
    ComputeMergedPolygons(const std::vector<Polygon>& polygons,
                          double minimalBuildingDistance)
    {
        // We merge the polygons by starting with the first vertex of each
        // polygon and walking counter-clockwise, adding either the next
        // vertex of the polygon or the polygon of another polygon if that
        // vertex is closer than a tolerance and generating a larger polygon.
        //
        // Note: This algorithm is O(n^2) and can be optimized by a more
        // clever search routine.

        // Make a copy of all polygons (so we can edit them in-place)
        std::vector<Polygon> mergedPolygons = polygons;

        // Create queue of polygons to check
        std::queue<size_t> polygonIndices;
        for (size_t i = 0; i < polygons.size(); i++)
            polygonIndices.push(i);

        // Check polygons until the queue is empty
        while (polygonIndices.size() > 0)
        {
            // Pop polygon from front of queue
            const size_t i = polygonIndices.front();
            polygonIndices.pop();

            std::cout << "Checking polygon " << i << std::endl;

            // Iterate over all other polygons
            for (size_t j = 0; j < polygons.size(); j++)
            {
                // Skip polygon itself
                if (i == j)
                    continue;

                // Skip if polygon has zero size (merged with other polygon)
                if (mergedPolygons[j].Points.size() == 0)
                    continue;

                // Compute squared distance between polygons
                const Polygon& Pi = mergedPolygons[i];
                const Polygon& Pj = mergedPolygons[j];
                const double d = Geometry::Distance2D(Pi, Pj);

                // Check if distance is smaller than the tolerance
                if (d < minimalBuildingDistance)
                {
                    std::cout << "CityModelGenerator: Buildings "
                              << i << " and " << j
                              << " are too close, merging." << std::endl;

                    // Compute merged polygon
                    Polygon mergedPolygon = MergedPolygons(Pi, Pj, d);

                    // Replace Pi, erase Pj and add Pi to queue
                    mergedPolygons[i] = mergedPolygon;
                    mergedPolygons[j].Points.clear();
                    polygonIndices.push(i);
                }
            }
        }

        return mergedPolygons;
    }

    static Polygon MergedPolygons(const Polygon& polygon0, const Polygon& p1,
                                  double distance)
    {
        // Compute intersection points for



    }

};

}

#endif
