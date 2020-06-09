// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_CITY_MODEL_GENERATOR_H
#define DTCC_CITY_MODEL_GENERATOR_H

#include <iostream>
#include <queue>
#include <vector>

#include "CityModel.h"
#include "GridField.h"
#include "Vector.h"
#include "Polygon.h"

namespace DTCC
{

class CityModelGenerator
{
public:
  // Generate city model from building footprints and height map
  static void GenerateCityModel(CityModel &cityModel,
                                const std::vector<Polygon> &polygons,
                                const GridField2D &heightMap,
                                double x0,
                                double y0,
                                double xMin,
                                double yMin,
                                double xMax,
                                double yMax)
  {
    std::cout << "CityModelGenerator: Generating city model..." << std::endl;

    // FIXME: Consider making all polygon processing in-place to
    // avoid copying the polygon data in each step.

    // Copy polygon data
    std::vector<Polygon> _polygons = polygons;

    // Compute transformed polygons
    _polygons =
        ComputeTransformedPolygons(_polygons, x0, y0, xMin, yMin, xMax, yMax);

    // Compute closed polygons
    _polygons = ComputeClosedPolygons(_polygons);

    // Compute counter-clockwise oriented polygons
    _polygons = ComputeOrientedPolygons(_polygons);

    // Add buildings
    for (auto const &polygon : _polygons)
    {
      Building building;
      building.Footprint = polygon;
      cityModel.Buildings.push_back(building);
    }

    // Compute building heights
    ComputeBuildingHeights(cityModel, heightMap);
  }

  // Simplify city model (simplify and merge polygons)
  static void SimplifyCityModel(CityModel &cityModel,
                                const GridField2D &heightMap,
                                double minimalBuildingDistance)
  {
    std::cout << "CityModelGenerator: Simplifying city model..." << std::endl;

    // FIXME: Consider making all polygon processing in-place to
    // avoid copying the polygon data in each step.

    // Copy polygon data
    std::vector<Polygon> _polygons;
    for (auto const &building : cityModel.Buildings)
      _polygons.push_back(building.Footprint);

    // Compute simplified polygons
    _polygons = ComputeSimplifiedPolygons(_polygons);

    // Compute merged polygons
    _polygons = ComputeMergedPolygons(_polygons, minimalBuildingDistance);

    // Compute simplified polygons (again)
    _polygons = ComputeSimplifiedPolygons(_polygons);

    // Add buildings
    cityModel.Buildings.clear();
    for (auto const &polygon : _polygons)
    {
      Building building;
      building.Footprint = polygon;
      cityModel.Buildings.push_back(building);
    }

    // Compute building heights (need to be recomputed for merged polygons)
    ComputeBuildingHeights(cityModel, heightMap);
  }

private:
  // Compute transformed polygons, keeping only polygons completely
  // within the domain and setting the origin to (xMin, yMin)
  static std::vector<Polygon>
  ComputeTransformedPolygons(const std::vector<Polygon> &polygons,
                             double x0,
                             double y0,
                             double xMin,
                             double yMin,
                             double xMax,
                             double yMax)
  {
    // Note that (x0, y0) are subtracted from the polygon coordinates.
    // The coordinates are then checked vs the domain dimensions.

    // Create empty list of polygons
    std::vector<Polygon> transformedPolygons;

    // Iterate over polygons
    for (auto const &polygon : polygons)
    {
      // Check if all points are inside
      bool inside = true;
      for (auto const &p : polygon.Vertices)
      {
        Vector2D q(p.x - x0, p.y - y0);
        if (q.x < xMin || q.y < yMin || q.x > xMax || q.y > yMax)
        {
          inside = false;
          break;
        }
      }

      // Add if inside
      if (inside)
      {
        Polygon transformedPolygon;
        for (auto const &p : polygon.Vertices)
        {
          Vector2D q(p.x - x0, p.y - y0);
          transformedPolygon.Vertices.push_back(q);
        }
        transformedPolygons.push_back(transformedPolygon);
      }
    }

    std::cout << "CityModelGenerator: Found " << transformedPolygons.size()
              << " building(s) out of " << polygons.size() << " inside domain"
              << std::endl;

    return transformedPolygons;
  }

  // Compute closed polygons (removing any duplicate vertices)
  static std::vector<Polygon>
  ComputeClosedPolygons(const std::vector<Polygon> &polygons)
  {
    // Create empty list of polygons
    std::vector<Polygon> closedPolygons;

    // Iterate over polygons
    for (auto const &polygon : polygons)
    {
      // Create empty polygon
      Polygon closedPolygon;

      // Iterate over points and add only unique points
      for (auto const &p : polygon.Vertices)
      {
        // Check if point is unique
        bool unique = true;
        for (auto const &q : closedPolygon.Vertices)
        {
          const double d = Geometry::Distance2D(p, q);
          if (d < Parameters::Epsilon)
          {
            unique = false;
            break;
          }
        }

        // Add point if unique
        if (unique)
          closedPolygon.Vertices.push_back(p);
      }

      // Add polygon
      closedPolygons.push_back(closedPolygon);
    }

    std::cout << "CityModelGenerator: Polygons closed" << std::endl;

    return closedPolygons;
  }

  // Compute counter-clockwise oriented polygons. It is assumed that the
  // input polygons are closed without duplicate vertices.
  static std::vector<Polygon>
  ComputeOrientedPolygons(const std::vector<Polygon> &polygons)
  {
    // Create empty list of polygons
    std::vector<Polygon> orientedPolygons;
    size_t numReversed = 0;

    // Iterate over polygons
    for (auto const &polygon : polygons)
    {
      // Copy polygon
      Polygon orientedPolygon = polygon;

      // Reverse polygon if not counter-clockwise
      if (Geometry::PolygonOrientation2D(orientedPolygon) != 0)
      {
        numReversed++;
        std::reverse(orientedPolygon.Vertices.begin(),
                     orientedPolygon.Vertices.end());
      }

      // Add polygon
      orientedPolygons.push_back(orientedPolygon);
    }

    std::cout << "CityModelGenerator: Reversed " << numReversed
              << " polygon(s) out of " << polygons.size() << std::endl;

    return orientedPolygons;
  }

  // Compute simplified polygons
  static std::vector<Polygon>
  ComputeSimplifiedPolygons(const std::vector<Polygon> &polygons)
  {
    // Create empty list of polygons
    std::vector<Polygon> simplifiedPolygons;

    // Iterate over polygons
    for (auto const &polygon : polygons)
    {
      // Previous vertex index
      const size_t numPoints = polygon.Vertices.size();
      size_t previousIndex = numPoints - 1;

      // Array of vertices to be included
      std::vector<bool> keepIndices(numPoints);
      std::fill(keepIndices.begin(), keepIndices.end(), true);

      // Create empty polygon
      Polygon simplifiedPolygon;

      // Iterate over vertices
      for (size_t i = 0; i < numPoints; i++)
      {
        // Get previous, current and next points
        const Vector2D &p0 = polygon.Vertices[previousIndex];
        const Vector2D &p1 = polygon.Vertices[i];
        const Vector2D &p2 = polygon.Vertices[(i + 1) % numPoints];

        // Compute angle (cosine)
        Vector2D u = p1 - p0;
        Vector2D v = p2 - p1;
        u /= u.Magnitude();
        v /= v.Magnitude();
        const double cos = Geometry::Dot2D(u, v);

        // Add vertex if angle is large enough
        if (cos < 1.0 - 0.1)
        {
          simplifiedPolygon.Vertices.push_back(p1);
          previousIndex = i;
        }
      }

      // Add simplified polygon
      simplifiedPolygons.push_back(simplifiedPolygon);
    }

    return simplifiedPolygons;
  }

  // Compute merged polygons. It is assumed that the input polygons
  // are closed without duplicate vertices and oriented.
  static std::vector<Polygon>
  ComputeMergedPolygons(const std::vector<Polygon> &polygons,
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
    for (size_t i = 0; i < mergedPolygons.size(); i++)
      polygonIndices.push(i);

    // Check polygons until the queue is empty
    while (polygonIndices.size() > 0)
    {
      // Pop polygon from front of queue
      const size_t i = polygonIndices.front();
      polygonIndices.pop();

      // Iterate over all other polygons
      for (size_t j = 0; j < mergedPolygons.size(); j++)
      {
        // Skip polygon itself
        if (i == j)
          continue;

        // Skip if polygon has zero size (merged with other polygon)
        if (mergedPolygons[j].Vertices.size() == 0)
          continue;

        // Compute squared distance between polygons
        const Polygon &Pi = mergedPolygons[i];
        const Polygon &Pj = mergedPolygons[j];
        const double d = Geometry::Distance2D(Pi, Pj);

        // Check if distance is smaller than the tolerance
        if (d < minimalBuildingDistance)
        {
          std::cout << "CityModelGenerator: Buildings " << i << " and " << j
                    << " are too close, merging" << std::endl;

          // Compute merged polygon
          Polygon mergedPolygon = MergePolygons(Pi, Pj, d);

          // Replace Pi, erase Pj and add Pi to queue
          mergedPolygons[i] = mergedPolygon;
          mergedPolygons[j].Vertices.clear();
          polygonIndices.push(i);
        }
      }
    }

    // Extract non-empty polygons
    std::vector<Polygon> _mergedPolygons;
    for (auto const &polygon : mergedPolygons)
    {
      if (polygon.Vertices.size() > 0)
        _mergedPolygons.push_back(polygon);
    }

    return _mergedPolygons;
  }

  // Merge the two polygons
  static Polygon MergePolygons(const Polygon &polygon0,
                               const Polygon &polygon1,
                               double distance)
  {
    // For now, we just compute the convex hull, consider
    // a more advanced merging later

    // Collect points
    std::vector<Vector2D> allPoints;
    for (auto const &p : polygon0.Vertices)
      allPoints.push_back(p);
    for (auto const &p : polygon1.Vertices)
      allPoints.push_back(p);

    // Remove duplicate points
    std::vector<Vector2D> uniquePoints;
    for (auto const &p : allPoints)
    {
      // Check if point is unique
      bool unique = true;
      for (auto const &q : uniquePoints)
      {
        const double d = Geometry::Distance2D(p, q);
        if (d < Parameters::Epsilon)
        {
          unique = false;
          break;
        }
      }

      // Add if unique
      if (unique)
        uniquePoints.push_back(p);
    }

    // Compute convex hull
    return Geometry::ConvexHull2D(uniquePoints);
  }

  // Compute building heights
  static void ComputeBuildingHeights(CityModel &cityModel,
                                     const GridField2D &heightMap)
  {
    // Iterate over buildings
    for (auto &building : cityModel.Buildings)
    {
      // Get building footprint
      Polygon &polygon = building.Footprint;

      // Compute center and radius of footprint
      const Vector2D center = Geometry::PolygonCenter2D(polygon);
      const double radius = Geometry::PolygonRadius2D(polygon, center);

      // Add points for sampling height
      std::vector<Vector2D> samplePoints;
      const double a = 0.5 * radius;
      samplePoints.push_back(center);
      samplePoints.push_back(Vector2D(center.x + a, center.y));
      samplePoints.push_back(Vector2D(center.x - a, center.y));
      samplePoints.push_back(Vector2D(center.x, center.y + a));
      samplePoints.push_back(Vector2D(center.x, center.y - a));

      // Compute mean height at points inside footprint
      double z = 0.0;
      size_t numInside = 0;
      for (auto const &p : samplePoints)
      {
        if (Geometry::PolygonContains2D(polygon, p))
        {
          z += heightMap(p);
          numInside += 1;
        }
      }

      // Check if we got at least one point
      if (numInside == 0)
      {
        std::cout << "CityModelGenerator: No sample points inside building, "
                     "setting height to 0"
                  << std::endl;
        numInside = 1;
      }

      // Set building height
      building.Height = z / numInside;
    }
  }
};

} // namespace DTCC

#endif
