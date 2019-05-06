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
                                  const HeightMap& heightMap)
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
        _polygons = ComputeMergedPolygons(_polygons);

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
    ComputeMergedPolygons(const std::vector<Polygon>& polygons)
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
            size_t i = polygonIndices.front();
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



            }

            // Pop polygon from queue

        }


        return mergedPolygons;

        // // Iterate over polygons
        // for (size_t m = 0; m < polygons.size(); m++)
        // {
        //     // Get current polygon
        //     const Polygon& Pm = polygons[m];

        //     // Create empty polygon
        //     Polygon mergedPolygon;

        //     // Add first vertex
        //     mergedPolygon.Points.push_back(Pm.Points[0]);

        //     // Iterate over remaining vertices
        //     for (size_t i = 0; i < Pm.Points.size(); i++)
        //     {
        //         // Get points and segment vector
        //         const Point2D& p0 = Pm.Points[i];
        //         const Point2D& p1 = Pm.Points[(i + 1) % Pm.Points.size()];
        //         const Point2D u = p1 - p0;

        //         // Iterate over all other polygons
        //         for (size_t n = 0; n < polygons.size(); n++)
        //         {
        //             // Skip the polygon itself
        //             if (m == n)
        //                 continue;

        //             // Get other polygon
        //             const Polygon& Pn = polygons[n];

        //             // Iterate over other polygon vertices
        //             for (size_t j = 0; j < Pn.Points.size(); j++)
        //             {
        //                 // Compute distance to segment

        //             }


        //         }


        //     }


    }

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
