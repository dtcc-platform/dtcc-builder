// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_CITY_MODEL_GENERATOR_H
#define DTCC_CITY_MODEL_GENERATOR_H

#include <iostream>
#include <queue>
#include <vector>

#include "CityModel.h"
#include "GridField.h"
#include "Plotting.h"
#include "PointCloud.h"
#include "Polyfix.h"
#include "Polygon.h"
#include "Timer.h"
#include "Vector.h"

namespace DTCC
{

class CityModelGenerator
{
public:
  /// Generate city model from building footprints, including only building
  /// inside the given bounding box. Note that this does not generated any
  /// building heights, only flat 2D buildings.
  ///
  /// @param cityModel The city model
  /// @param footprints Footprints of buildings (polygons)
  /// @param UUIDs UUIDs of buildings
  /// @param entityIDs Indices of buildings (in shapefile)
  /// @param bbox Bounding box of domain
  static void GenerateCityModel(CityModel &cityModel,
                                const std::vector<Polygon> &footprints,
                                const std::vector<std::string> &UUIDs,
                                const std::vector<int> &entityIDs,
                                const BoundingBox2D &bbox)
  {
    Info("CityModelGenerator: Generating city model...");
    Timer("GenerateCityModel");

    // Clear data
    cityModel.Buildings.clear();

    // Add buildings
    for (size_t i = 0; i < footprints.size(); i++)
    {
      // Skip if not inside  bounding box
      if (!Geometry::BoundingBoxContains2D(bbox, footprints[i]))
        continue;

      // Add building
      Building building;
      building.Footprint = footprints[i];
      building.UUID = UUIDs[i];
      building.SHPFileID = entityIDs[i];
      cityModel.Buildings.push_back(building);
    }

    Info("CityModelGenerator: Added " + str(cityModel.Buildings.size()) + "/" +
         str(footprints.size()) + " inside bounding box");
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

  /// Extract roof points for buildings from a point cloud.
  /// The roof points of a building are defined as all points
  /// that fall within the footprint of a building. The points
  /// are ordered by increasing height (z-coordinate).
  ///
  /// @param cityModel The city model
  /// @param pointCloud Point cloud
  static void ExtractRoofPoints(CityModel &cityModel,
                                const PointCloud &pointCloud)
  {
    Info("CityModelGenerator: Extracting roof points...");
    Timer("ExtractRoofPoints");

    // Build search trees
    pointCloud.BuildSearchTree();
    cityModel.BuildSearchTree();

    // Compute bounding box tree collisions
    const auto collisions = pointCloud.bbtree.Find(cityModel.bbtree);

    // Clear previous roof points (if any)
    for (auto &building : cityModel.Buildings)
      building.RoofPoints.clear();

    // Add points to buildings. Note that we only which points
    // belong to the bounding boxes of the buildings so we need
    // to check for each point whether it belongs to the footprint
    // of the building.
    for (auto &index : collisions)
    {
      // Get point and building
      const Point3D &point = pointCloud.Points[index.first];
      Building &building = cityModel.Buildings[index.second];

      // Check if point is inside building footprint
      const Point2D p2D{point.x, point.y};
      if (Geometry::PolygonContains2D(building.Footprint, p2D))
        building.RoofPoints.push_back(point);
    }

    // Sort building points by height
    for (auto &building : cityModel.Buildings)
    {
      std::vector<Point3D> &points{building.RoofPoints};
      std::sort(
          points.begin(), points.end(),
          [](const Point3D &p, const Point3D &q) -> bool { return p.z < q.z; });
    }

    // Uncomment for debugging (plot footprints and points)
    /*
    Plotting::Init();
    for (size_t i = 0; i < cityModel.Buildings.size(); i++)
    {
      const auto building = cityModel.Buildings[i];
      const auto points = building.RoofPoints;
      Info("Building " + str(i) + ": n = " + str(points.size()));
      if (points.size() > 0)
      {
        Plotting::Plot(building);
        Plotting::Plot(points);
      }
    }
    */

    // Compute some statistics
    size_t min = std::numeric_limits<size_t>::max();
    size_t max = 0;
    size_t sum = 0;
    for (auto &building : cityModel.Buildings)
    {
      const size_t n = building.RoofPoints.size();
      min = std::min(min, n);
      max = std::max(max, n);
      sum += n;
    }
    const double mean = static_cast<double>(sum) /
                        static_cast<double>(cityModel.Buildings.size());

    Info("CityModelGenerator: min/mean/max number of points per building is " +
         str(min) + "/" + str(mean) + "/" + str(max));
  }

  /// Compute heights of buildings from roof points. This requires
  /// that ExtractRoofPoints() has been called to extract roof points
  /// from point cloud data.
  ///
  /// This version produces building heights of better quality since
  /// it avoids smoothing effect that is introduced when computing
  /// heights from the DSM.
  ///
  /// @param cityModel The city model
  /// @param dtm Digital Terrain Model (excluding buildings)
  /// @param heightPercentile Percentile to use for setting building height
  static void ComputeBuildingHeights(CityModel &cityModel,
                                     const GridField2D &dtm,
                                     double heightPercentile)
  {
    Info("CityModelGenerator: Computing building heights...");
    Timer("ComputeBuildingHeights");

    // Iterate over buildings
    for (auto &building : cityModel.Buildings)
    {
      // Skip if missing points
      std::vector<Point3D> &points{building.RoofPoints};
      if (points.empty())
      {
        Warning("Missing roof points for building " + building.UUID);
        continue;
      }

      // Compute absolute height of building from roof points
      size_t index = std::max(0.0, heightPercentile * points.size());
      index = std::min(points.size(), index);
      const double h0 = points[index].z;

      // Compute absolute height of ground from DTM (sample at center of
      // building)
      const Point2D center = Geometry::PolygonCenter2D(building.Footprint);
      const double h1 = dtm(center);

      // Set building height(s)
      building.Height = h1 - h0;
      building.GroundHeight = h0;
    }
  }

  /// Compute heights of buildings from DSM data.
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

  /// Generate a random city model. Used for benchmarking.
  ///
  /// @param cityModel The city model
  /// @param dtm Digital Terrian Map
  /// @param numBuildings Number of buildings
  static void RandomizeCityModel(CityModel &cityModel,
                                 const GridField2D &dtm,
                                 size_t numBuildings)
  {
    Info("CityModelGenerator: Randomizing city model...");

    // Some hard-coded building dimensions
    const double A = 20.0;  // Maximum building side length
    const double H = 10.0;  // Maximum building height
    const double N = 10000; // Maximum number of attempts

    // Get bounding box of domain
    const BoundingBox2D &bbox = dtm.Grid.BoundingBox;
    const double dx = bbox.Q.x - bbox.P.x;
    const double dy = bbox.Q.y - bbox.P.y;

    // Iterate over the number of buildings to generate
    std::vector<Point2D> centers;
    for (size_t i = 0; i < numBuildings; i++)
    {
      // Keep trying until we find an empty spot
      size_t counter = 0;
      while (true)
      {
        // Check number of attempts
        if (++counter > N)
        {
          Info("Try setting a smaller number of random buildings.");
          Error("Unable to randomize city model; reached maximum number of "
                "attempts.");
        }

        // Randomize center of building
        Point2D c(bbox.P.x + Utils::Random() * dx,
                  bbox.P.y + Utils::Random() * dy);

        // Check that we are not too close to other buildings, but
        // note that buildings may actually overlap slightly which may
        // also happen with real-world data.
        bool ok = true;
        for (auto const &p : centers)
        {
          const double d = Geometry::Distance2D(p, c);
          if (d < 0.5 * A)
          {
            ok = false;
            break;
          }
        }
        if (!ok)
          continue;

        // Check that we are not close to the domain boundary
        if ((c.x - bbox.P.x < 2.0 * A) || (bbox.Q.x - c.x < 2.0 * A) ||
            (c.y - bbox.P.y < 2.0 * A) || (bbox.Q.y - c.y < 2.0 * A))
          continue;

        // Randomize dimension
        const double a = (0.05 + 0.95 * Utils::Random()) * A;
        const double b = (0.05 + 0.95 * Utils::Random()) * A;
        const double h = (0.25 + 0.75 * Utils::Random()) * H;

        // Generate building
        Building building = GenerateBuilding(c, a, b, h, dtm(c));

        // Add building
        cityModel.Buildings.push_back(building);
        centers.push_back(c);

        Info("Creating random building " + str(i + 1) + "/" +
             str(numBuildings) + " at c = " + str(c));
        break;
      }
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

  // Generate building with given dimensions
  static Building GenerateBuilding(
      const Point2D &c, double a, double b, double height, double groundHeight)
  {
    Building building;
    building.Footprint.Vertices.push_back(
        Point2D{c.x - 0.5 * a, c.y - 0.5 * b});
    building.Footprint.Vertices.push_back(
        Point2D{c.x + 0.5 * a, c.y - 0.5 * b});
    building.Footprint.Vertices.push_back(
        Point2D{c.x + 0.5 * a, c.y + 0.5 * b});
    building.Footprint.Vertices.push_back(
        Point2D{c.x - 0.5 * a, c.y + 0.5 * b});
    building.Height = height;
    building.GroundHeight = groundHeight;
    return building;
  }
};

} // namespace DTCC

#endif
