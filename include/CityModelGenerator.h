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
#include "PointCloudProcessor.h"
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
         str(footprints.size()) + " buildings inside bounding box");
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

  /// Extract ground and roof points from point cloud.
  ///
  /// The ground points of a building are defined as all points
  /// of class 2 (Ground) or 9 (Water) that fall within a given
  /// distance from the building footprint.
  ///
  /// The roof points of a building are defined as all points
  /// of class 6 (Building) that fall within the building
  /// footprint. However, since that classification seems to
  /// be missing in the data from LM, we are currently using
  /// all points (except class 2 and 9).
  ///
  /// @param cityModel The city model
  /// @param pointCloud Point cloud (unfiltered)
  /// @param groundMargin Margin around building for detecting ground points
  static void ExtractBuildingPoints(CityModel &cityModel,
                                    const PointCloud &pointCloud,
                                    double groundMargin)
  {
    Info("CityModelGenerator: Extracting building points...");
    Timer("ExtractBuildingPoints");

    // Check that point cloud is not empty
    if (pointCloud.Points.empty())
      Error("Empty point cloud");

    // Check that point cloud has classifications
    if (pointCloud.Points.size() != pointCloud.Classification.size())
      Error("Missing classifications for point cloud");

    // Build search trees
    pointCloud.BuildSearchTree(true);
    cityModel.BuildSearchTree(true, groundMargin);

    // Compute bounding box tree collisions
    const auto collisions = pointCloud.bbtree.Find(cityModel.bbtree);

    // Clear old building points (if any)
    for (auto &building : cityModel.Buildings)
    {
      building.GroundPoints.clear();
      building.RoofPoints.clear();
    }

    // Squared margin for detecting ground points
    const double d2 = groundMargin * groundMargin;

    // Iterate over collisions and extract points
    for (auto &index : collisions)
    {
      // Get point and building
      const Point3D &p3D = pointCloud.Points[index.first];
      const Point2D p2D{p3D.x, p3D.y};
      const uint8_t clf = pointCloud.Classification[index.first];
      Building &building = cityModel.Buildings[index.second];

      // Check for ground points
      if (clf == 2 || clf == 9)
      {
        if (Geometry::SquaredDistance2D(building.Footprint, p2D) < d2)
          building.GroundPoints.push_back(p3D);
      }

      // Check for roof points
      // else if (clf == 6)
      else
      {
        if (Geometry::PolygonContains2D(building.Footprint, p2D))
          building.RoofPoints.push_back(p3D);
      }
    }

    // Sort points by height
    for (auto &building : cityModel.Buildings)
    {
      std::sort(
          building.GroundPoints.begin(), building.GroundPoints.end(),
          [](const Point3D &p, const Point3D &q) -> bool { return p.z < q.z; });
      std::sort(
          building.RoofPoints.begin(), building.RoofPoints.end(),
          [](const Point3D &p, const Point3D &q) -> bool { return p.z < q.z; });
    }

    // Compute some statistics
    size_t minG{std::numeric_limits<size_t>::max()};
    size_t minR{std::numeric_limits<size_t>::max()};
    size_t maxG{0}, maxR{0};
    size_t sumG{0}, sumR{0};
    for (const auto &building : cityModel.Buildings)
    {
      // Ground points
      const size_t nG = building.GroundPoints.size();
      minG = std::min(minG, nG);
      maxG = std::max(maxG, nG);
      sumG += nG;

      // Roof points
      const size_t nR = building.RoofPoints.size();
      minR = std::min(minR, nR);
      maxR = std::max(maxR, nR);
      sumR += nR;
    }
    const double meanG = static_cast<double>(sumG) / cityModel.Buildings.size();
    const double meanR = static_cast<double>(sumR) / cityModel.Buildings.size();

    Info("CityModelGenerator: min/mean/max number of ground points per "
         "building is " +
         str(minG) + "/" + str(meanG) + "/" + str(maxG));
    Info("CityModelGenerator: min/mean/max number of roof points per building "
         "is " +
         str(minR) + "/" + str(meanR) + "/" + str(maxR));
  }

  /// Compute heights of buildings from ground and roof points. This
  /// requires that ExtractBuildingPoints() has been called to extract
  /// the points from point cloud data.
  ///
  /// @param cityModel The city model
  /// @param groundPercentile Percentile used for setting ground height
  /// @param roofPercentile Percentile used for setting roof height
  static void ComputeBuildingHeights(CityModel &cityModel,
                                     double groundPercentile,
                                     double roofPercentile)
  {
    Info("CityModelGenerator: Computing building heights...");
    Timer("ComputeBuildingHeights");

    Plotting::Init();

    // Count the number of successful buildings
    size_t numSuccess = 0;

    // Iterate over buildings
    for (auto &building : cityModel.Buildings)
    {
      // Skip if missing points
      if (building.GroundPoints.empty() || building.RoofPoints.empty())
      {
        Warning("Missing points for building " + building.UUID);
        continue;
      }

      // Compute heights
      double h0 = GetPercentile(building.GroundPoints, groundPercentile).z;
      double h1 = GetPercentile(building.RoofPoints, roofPercentile).z;

      // Check that h0 < h1
      if (!(h0 < h1))
      {
        Warning("Roof height lower than ground height for building " +
                building.UUID);
        h1 = h0;
        continue;
      }

      // Set building height(s)
      building.Height = h1 - h0;
      building.GroundHeight = h0;
      numSuccess++;

      // Uncomment for debugging
      Plotting::Plot(building);
    }

    // Report number of failed buildings
    const size_t numFail = cityModel.Buildings.size() - numSuccess;
    Info("CityModelGenerator: Height computation failed for " + str(numFail) +
         "/" + str(cityModel.Buildings.size()) + " buildings");
  }

  // Generate a random city model. Used for benchmarking.
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
  // Get percentile object from array. It is assumed that the array is ordered.
  template <class T>
  static T GetPercentile(const std::vector<T> &array, double percentile)
  {
    size_t index = std::max(0.0, percentile * array.size());
    index = std::min(index, array.size() - 1);
    return array[index];
  }

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
