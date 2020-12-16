// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_CITY_MODEL_GENERATOR_H
#define DTCC_CITY_MODEL_GENERATOR_H

#include <iostream>
#include <queue>
#include <vector>

#include "CityModel.h"
#include "GridField.h"
#include "Polyfix.h"
#include "Polygon.h"
#include "Timer.h"
#include "Vector.h"

namespace DTCC
{

class CityModelGenerator
{
public:
  /// Generate city model from building footprints. The given origin is
  /// subtracted from all coordinates and only buildings completely inside
  /// the given bounding box (after coordinate transformation) are included.
  ///
  /// @param cityModel The city model
  /// @param footprints Building footprints (polygons)
  /// @param origin Origin to be subtracted
  /// @param bbox Bounding box of domain
  static void GenerateCityModel(CityModel &cityModel,
                                const std::vector<Polygon> &footprints,
                                const std::vector<std::string> &UUIDs,
                                const std::vector<int> &entityIDs,
                                const Point2D &origin,
                                const BoundingBox2D &bbox)
  {
    Info("CityModelGenerator: Generating city model...");
    Timer("GenerateCityModel");

    // Clear old data
    cityModel.Buildings.clear();

    // Add buildings
    // for (const auto &footprint : footprints)
    for (size_t i = 0; i < footprints.size(); i++)
    {
      // Create transformed footprint
      Polygon transformedFootprint = footprints[i];
      Polyfix::Transform(transformedFootprint, origin);

      // Add if inside bounding box
      if (Geometry::BoundingBoxContains2D(bbox, transformedFootprint))
      {
        Building building;
        building.Footprint = transformedFootprint;
        building.UUID = UUIDs[i];
        // Uncomment for debugging
        // building.debugID = i + 1;
        // Add SHP file entityID
        building.SHPFileID = entityIDs[i];
        cityModel.Buildings.push_back(building);
        // std::cout << "i = " << i << " entityID = " << entityIDs[i] << " UUID
        // = " << UUIDs[i] << std::endl;
      }
    }

    Info("CityModelGenerator: Found " + str(cityModel.Buildings.size()) + "/" +
         str(footprints.size()) + " buildings inside domain");
  }

  /// Clean city model by making sure that all building footprints
  /// are closed and counter-clockwise oriented.
  ///
  /// @param cityModel The city model
  /// @param minimalVertexDistance Minimal vertex distance
  static void CleanCityModel(CityModel &cityModel, double minimalVertexDistance)
  {
    Info("CityModelGenerator: Cleaning city model...");
    Timer("CleanCityModel");

    // Iterate over buildings
    size_t numClosed = 0;
    size_t numOriented = 0;
    size_t numVertexMerged = 0;
    size_t numEdgeMerged = 0;
    for (auto &building : cityModel.Buildings)
    {
      // Make closed
      numClosed += Polyfix::MakeClosed(building.Footprint, Parameters::Epsilon);

      // Make oriented
      numOriented += Polyfix::MakeOriented(building.Footprint);

      // Merge vertices
      numVertexMerged +=
          Polyfix::MergeVertices(building.Footprint, minimalVertexDistance);

      // Merge edges
      numEdgeMerged += Polyfix::MergeEdges(building.Footprint,
                                           Parameters::FootprintAngleThreshold);
    }

    Info("CityModelGenerator: Fixed " + str(numClosed) + "/" +
         str(cityModel.Buildings.size()) + " polygons that were not closed");
    Info("CityModelGenerator: Fixed " + str(numOriented) + "/" +
         str(cityModel.Buildings.size()) + " polygons that were not oriented");
    Info("CityModelGenerator: Merged vertices for " + str(numVertexMerged) +
         "/" + str(cityModel.Buildings.size()) + " polygons");
    Info("CityModelGenerator: Merged edges for " + str(numEdgeMerged) + "/" +
         str(cityModel.Buildings.size()) + " polygons");
  }

  // Simplify city model (simplify and merge polygons)
  static void SimplifyCityModel(CityModel &cityModel,
                                double minimalBuildingDistance,
                                double minimalVertexDistance)
  {
    Info("CityModelGenerator: Simplifying city model...");
    Timer("SimplifyCityModel");

    // Merge buildings if too close
    MergeBuildings(cityModel, minimalBuildingDistance);

    // Clean after merge
    CleanCityModel(cityModel, minimalVertexDistance);
  }

  /// Compute heights of buildings from height map.
  ///
  /// @param cityModel The city model
  /// @param dsm       Digital Surface Model (including buildings)
  /// @param dtm       Digital Terrain Model (excluding buildings)
  static void ComputeBuildingHeights(CityModel &cityModel,
                                     const GridField2D &dsm,
                                     const GridField2D &dtm)
  {
    Info("CityModelGenerator: Computing building heights...");
    Timer("ComputeBuildingHeights");

    // Iterate over buildings
    for (size_t i = 0; i < cityModel.Buildings.size(); i++)
    {
      // Get building
      Building &building = cityModel.Buildings[i];

      // Compute center and radius of building footprint
      const Point2D center = Geometry::PolygonCenter2D(building.Footprint);
      const double radius =
          Geometry::PolygonRadius2D(building.Footprint, center);

      // Add points for sampling height
      std::vector<Point2D> samplePoints;
      samplePoints.push_back(center);
      const size_t m = 2;
      const size_t n = 8;
      for (size_t k = 1; k <= m; k++)
      {
        const double r =
            static_cast<double>(k) / static_cast<double>(m) * 0.5 * radius;
        for (size_t l = 0; l < n; l++)
        {
          const double theta =
              static_cast<double>(l) / static_cast<double>(n) * 2.0 * M_PI;
          Point2D p(center.x + r * cos(theta), center.y + r * sin(theta));
          samplePoints.push_back(p);
        }
      }

      // Compute mean height at points inside footprint
      double h0 = 0.0;
      double h1 = 0.0;
      size_t numInside = 0;
      for (auto const &p : samplePoints)
      {
        if (Geometry::PolygonContains2D(building.Footprint, p))
        {
          h0 += dtm(p);
          h1 += dsm(p);
          numInside += 1;
        }
      }

      // Check if we got at least one point
      if (numInside == 0)
      {
        Info("CityModelGenerator: No sample points inside building " + str(i) +
             ", setting height to 0");
        numInside = 1;
      }

      // Compute mean
      h0 /= static_cast<double>(numInside);
      h1 /= static_cast<double>(numInside);

      // Set building height(s)
      building.Height = h1 - h0;
      building.GroundHeight = h0;
    }
  }

private:
  // Merge all buildings closer than a given distance
  static void MergeBuildings(CityModel &cityModel,
                             double minimalBuildingDistance)
  {
    Info("CityModelGenerator: Merging buildings...");

    // Avoid using sqrt for efficiency
    const double tol2 = minimalBuildingDistance * minimalBuildingDistance;

    // Get buildings
    std::vector<Building> &buildings = cityModel.Buildings;

    // Create queue of indices to check
    std::queue<size_t> indices;
    for (size_t i = 0; i < buildings.size(); i++)
      indices.push(i);

    // Process queue until empty
    size_t numMerged = 0;
    while (!indices.empty())
    {
      // Pop index of next building to check
      const size_t i = indices.front();
      indices.pop();

      // Iterate over all other buildings
      for (size_t j = 0; j < buildings.size(); j++)
      {
        // Skip building itself
        if (i == j)
          continue;

        // Skip if merged with other building (size set to 0)
        if (buildings[j].Footprint.Vertices.empty())
          continue;

        // Compute squared distance between polygons
        const Polygon &Pi = buildings[i].Footprint;
        const Polygon &Pj = buildings[j].Footprint;
        const double d2 = Geometry::SquaredDistance2D(Pi, Pj);

        // Check if distance is smaller than the tolerance
        if (d2 < tol2)
        {
          Progress("CityModelGenerator: Buildings " + str(i) + " and " +
                   str(j) + " are too close, merging");

          // Compute merged polygon
          Polygon mergedPolygon =
              Polyfix::MergePolygons(Pi, Pj, minimalBuildingDistance);
          numMerged++;

          // Replace Pi, erase Pj and add Pi to queue
          buildings[i].Footprint = mergedPolygon;
          buildings[j].Footprint.Vertices.clear();
          indices.push(i);
        }
      }
    }

    // Extract non-empty polygons
    std::vector<Building> mergedBuildings;
    for (const auto &building : buildings)
    {
      if (!building.Footprint.Vertices.empty())
        mergedBuildings.push_back(building);
    }

    // Overwrite buildings
    cityModel.Buildings = mergedBuildings;

    Info("CityModelGenerator: Merged " + str(numMerged) + " buildings");
  }

  // Merge the two polygons
  static Polygon MergePolygons(const Polygon &polygon0, const Polygon &polygon1)
  {
    // For now, we just compute the convex hull, consider
    // a more advanced merging later

    // Collect points
    std::vector<Point2D> allPoints;
    for (auto const &p : polygon0.Vertices)
      allPoints.push_back(p);
    for (auto const &p : polygon1.Vertices)
      allPoints.push_back(p);

    // Remove duplicate points
    std::vector<Point2D> uniquePoints;
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
};

} // namespace DTCC

#endif
