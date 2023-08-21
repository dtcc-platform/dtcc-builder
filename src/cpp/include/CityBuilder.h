// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_CITY_BUILDER_H
#define DTCC_CITY_BUILDER_H

#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <unordered_set>
#include <vector>

#include "BuildingProcessor.h"
#include "KDTreeVectorOfVectorsAdaptor.h"
#include "Logging.h"
#include "PointCloudProcessor.h"
#include "Polyfix.h"
#include "Timer.h"
#include "model/City.h"
#include "model/GridField.h"
#include "model/PointCloud.h"
#include "model/Polygon.h"
#include "model/Vector.h"

namespace DTCC_BUILDER
{

class CityBuilder
{
public:
  /// Build city from building footprints, including only building
  /// inside the given bounding box. Note that this does not generate any
  /// building heights, only flat 2D buildings.
  ///
  /// @param city The city
  /// @param footprints Footprints of buildings (polygons)
  /// @param UUIDs UUIDs of buildings
  /// @param entityIDs Indices of buildings (in shapefile)
  /// @param bbox Bounding box of domain
  /// @param minBuildingDistance Minimal distance from building to domain
  /// boundary
  static void BuildCity(City &city,
                        const std::vector<Polygon> &footprints,
                        const std::vector<std::string> &UUIDs,
                        const std::vector<int> &entityIDs,
                        const BoundingBox2D &bbox,
                        double minBuildingDistance,
                        double minBuildingSize)
  {
    info("Building city...");
    Timer timer("BuildCity");

    // Clear data
    city.Buildings.clear();

    // Add buildings
    for (size_t i = 0; i < footprints.size(); i++)
    {
      // Skip if not inside  bounding box
      if (!Geometry::BoundingBoxContains2D(bbox, footprints[i],
                                           minBuildingDistance))
      {
        warning("Skipping building " + UUIDs[i] +
                "; outside domain or too close to boundary");
        continue;
      }

      // Add building
      Building building;
      building.Footprint = footprints[i];
      building.UUID = UUIDs[i];
      building.SHPFileID = entityIDs[i];
      if (Geometry::PolygonArea(building.Footprint) < minBuildingSize)
      {
        building.error |= BuildingError::BUILDING_TOO_SMALL;
      }

      city.Buildings.push_back(building);
    }

    info("Added " + str(city.Buildings.size()) + "/" + str(footprints.size()) +
         " buildings inside bounding box");
  }

  /// Clean city by making sure that all building footprints
  /// are closed and counter-clockwise oriented.
  ///
  /// @param city The city
  /// @param min_vertex_distance Minimal vertex distance
  static City clean_city(const City &city, double min_vertex_distance)
  {
    info("Cleaning city...");
    Timer timer("clean_city");

    // Clear search tree (since it might become invalid)
    city.bbtree.Clear();

    // Count some stats
    size_t numClosed = 0;
    size_t numOriented = 0;
    size_t numVertexMerged = 0;
    size_t numEdgeMerged = 0;
    size_t numRemoved = 0;

    // Initialize new city
    City _city;
    _city.Name = city.Name;

    // Clean buildings
    for (auto &building : city.Buildings)
    {
      // Make copy of building
      Building _building{building};

      // Make closed
      numClosed += Polyfix::MakeClosed(_building.Footprint, Constants::Epsilon);

      // Make oriented
      numOriented += Polyfix::MakeOriented(_building.Footprint);

      // Merge vertices (but skip if only 4 vertices or less)
      if (_building.Footprint.Vertices.size() > 4)
      {
        numVertexMerged +=
            Polyfix::MergeVertices(_building.Footprint, min_vertex_distance);
      }

      // Merge edges (but skip if only 4 vertices or less)
      if (_building.Footprint.Vertices.size() > 4)
      {
        numEdgeMerged += Polyfix::MergeEdges(
            _building.Footprint, Constants::FootprintAngleThreshold);
      }

      // Keep only valid buildings
      if (_building.Valid())
        _city.Buildings.push_back(_building);
      else
        numRemoved++;
    }

    info("Fixed " + str(numClosed) + "/" + str(city.Buildings.size()) +
         " polygons that were not closed");
    info("Fixed " + str(numOriented) + "/" + str(city.Buildings.size()) +
         " polygons that were not oriented");
    info("Merged vertices for " + str(numVertexMerged) + "/" +
         str(city.Buildings.size()) + " polygons");
    info("Merged edges for " + str(numEdgeMerged) + "/" +
         str(city.Buildings.size()) + " polygons");
    info("Removed " + str(numRemoved) + "/" + str(city.Buildings.size()) +
         " buildings (invalid/too small after cleaning)");

    return _city;
  }

  /// Compute ground and roof points from point cloud.
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
  /// @param city The city
  /// @param pointCloud Point cloud (unfiltered)
  /// @param groundMargin Margin around building for detecting ground points
  static City compute_building_points(const City &city,
                                      const PointCloud &pointCloud,
                                      double groundMargin,
                                      double groundOutlierMargin)

  {
    info("Computing building points...");
    Timer timer("compute_building_points");

    // Check that point cloud is not empty
    if (pointCloud.Points.empty())
      error("Empty point cloud");

    // Check that point cloud has classifications
    if (pointCloud.Points.size() != pointCloud.Classifications.size())
      error("Missing classifications for point cloud");

    bool classifiedPoints =
        (pointCloud.Points.size() == pointCloud.Classifications.size());
    bool classifedBuildings = false;
    if (classifiedPoints)
      classifedBuildings = pointCloud.HasClassification(6);

    auto kdt_timer = Timer("ExtractBuildingPoints: BuildKDTree");
    // build a kd-tree for radius search
    typedef KDTreeVectorOfVectorsAdaptor<std::vector<Point3D>, double,
                                         2 /* dims */>
        my_kd_tree_t;
    my_kd_tree_t pc_index(2, pointCloud.Points, 20 /* max leaf */);
    kdt_timer.Stop();

    // Create copy of city
    City _city{city};

    // Iterate over buildings
    for (auto &building : _city.Buildings)
    {
      building.GroundPoints.clear();
      building.RoofPoints.clear();

      auto centerPoint = Geometry::PolygonCenter2D(building.Footprint);
      double radius =
          Geometry::PolygonRadius2D(building.Footprint, centerPoint);
      radius *= radius;
      radius += groundMargin;

      std::vector<double> query_pt{centerPoint.x, centerPoint.y};
      auto radius_t = Timer("RadiusQuery");
      auto indices_dists = pc_index.radiusQuery(&query_pt[0], radius);
      radius_t.Stop();
      for (auto const &ind_pt : indices_dists)
      {
        size_t idx = ind_pt.first;
        const uint8_t clf = pointCloud.Classifications[idx];
        const Point3D &p3D = pointCloud.Points[idx];
        const Point2D p2D{p3D.x, p3D.y};

        if (classifiedPoints)
        {
          if (clf == 2 || clf == 9)
          {
            building.GroundPoints.push_back(p3D);
          }
          else if (clf == 6 || (!classifedBuildings && clf < 2))
          {
            // auto pc_timer = Timer("PolygoCOntains2D");
            if (Geometry::PolygonContains2D(building.Footprint, p2D))
            {
              building.RoofPoints.push_back(p3D);
            }
            // pc_timer.Stop();
          }
        }
        else // unclassified data
        {
          if (Geometry::PolygonContains2D(building.Footprint, p2D))
          {
            building.RoofPoints.push_back(p3D);
          }
          else
          {
            // all points not in the roof polygon are considered ground
            building.GroundPoints.push_back(p3D);
          }
        }
      }
    }

    // Remove ground outliers
    size_t numGroundPoints = 0;
    size_t numGroundOutliers = 0;
    for (auto &building : _city.Buildings)
    {
      // Count total number of points
      numGroundPoints += building.GroundPoints.size();

      // Remove outliers and count total number of outliers
      numGroundOutliers += PointCloudProcessor::RemoveOutliers(
                               building.GroundPoints, groundOutlierMargin)
                               .size();
    }
    const double outlierGroundPercentage =
        (100.0 * numGroundOutliers) / numGroundPoints;
    info("Removed ground point outliers (" + str(outlierGroundPercentage) +
         "%)");

    double ptsPrSqm;
    double pointCoverage;
    size_t tooFew = 0;
    for (auto &building : _city.Buildings)
    {
      ptsPrSqm = static_cast<double>(building.RoofPoints.size()) /
                 Geometry::PolygonArea(building.Footprint);
      if (ptsPrSqm < 0.25)
      {
        building.error |= BuildingError::BUILDING_TOO_FEW_POINTS;
        tooFew++;
      }
      pointCoverage = BuildingProcessor::PointCoverage(building, 2.0);
      // info("PointCoverage: " + str(pointCoverage));
      if (pointCoverage < 0.5)
      {
        building.error |= BuildingError::BUILDING_INSUFFICIENT_POINT_COVERAGE;
      }
    }
    info("Number of buildings with too few roof points: " + str(tooFew));

    // Sort points by height
    for (auto &building : _city.Buildings)
    {
      std::sort(building.GroundPoints.begin(), building.GroundPoints.end(),
                [](const Point3D &p, const Point3D &q) -> bool
                { return p.z < q.z; });
      std::sort(building.RoofPoints.begin(), building.RoofPoints.end(),
                [](const Point3D &p, const Point3D &q) -> bool
                { return p.z < q.z; });
    }

    // Compute some statistics
    size_t minG{std::numeric_limits<size_t>::max()};
    size_t minR{std::numeric_limits<size_t>::max()};
    size_t maxG{0}, maxR{0};
    size_t sumG{0}, sumR{0};
    for (const auto &building : _city.Buildings)
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
    const double meanG = static_cast<double>(sumG) / _city.Buildings.size();
    const double meanR = static_cast<double>(sumR) / _city.Buildings.size();

    info("min/mean/max number of ground points per "
         "building is " +
         str(minG) + "/" + str(meanG) + "/" + str(maxG));
    info("min/mean/max number of roof points per building "
         "is " +
         str(minR) + "/" + str(meanR) + "/" + str(maxR));

    return _city;
  }

  /// Compute heights of buildings from ground and roof points. This
  /// requires that ExtractBuildingPoints() has been called to extract
  /// the points from point cloud data.
  ///
  /// @param city The city
  /// @param dtm Digital Terrain Map, used for ground height if points are
  /// missing
  /// @param groundPercentile Percentile used for setting ground height
  /// @param roofPercentile Percentile used for setting roof height
  static void ComputeBuildingHeights(City &city,
                                     const GridField &dtm,
                                     double groundPercentile,
                                     double roofPercentile)
  {
    info("Computing building heights...");
    Timer timer("ComputeBuildingHeights");

    // FIXME: Make this a parameter?
    // FIXME: How do we treat this in relation to layer height?
    // FIXME: Also a problem inside MeshGenerator (warning in TrimMesh3D)
    // FIXME: Make this a parameter
    const double minBuildingHeight{2.5};

    // Uncomment for debugging
    // Plotting::Init();

    // Count missing or bad data
    size_t numMissingGroundPoints = 0;
    size_t numMissingRoofPoints = 0;
    size_t numSmallHeights = 0;

    // Iterate over buildings
    for (auto &building : city.Buildings)
    {
      // Compute ground height h0
      double h0{0};
      if (building.GroundPoints.empty())
      {
        // warning("Missing ground points for building " + building.UUID);
        // info("Setting ground height from DTM");
        h0 = dtm(Geometry::PolygonCenter2D(building.Footprint));
        numMissingGroundPoints++;
        building.error |= BuildingError::BUILDING_NO_GROUND_POINTS;
      }
      else
      {
        // Pick percentile from ground points
        sort(building.GroundPoints.begin(), building.GroundPoints.end(),
             [](const Point3D &lhs, const Point3D &rhs)
             { return lhs.z < rhs.z; });
        h0 = GetPercentile(building.GroundPoints, groundPercentile).z;
      }

      // Compute roof height h1
      double h1{0};
      if (building.RoofPoints.empty())
      {
        // warning("Missing roof points for building " + building.UUID);
        // info("Setting building height to " +
        //     str(minBuildingHeight) + "m");
        h1 = h0 + minBuildingHeight;
        numMissingRoofPoints++;
        building.error |= BuildingError::BUILDING_NO_ROOF_POINTS;
      }
      else
      {
        sort(building.RoofPoints.begin(), building.RoofPoints.end(),
             [](const Point3D &lhs, const Point3D &rhs)
             { return lhs.z < rhs.z; });
        h1 = GetPercentile(building.RoofPoints, roofPercentile).z;
      }

      // Check that h0 < h1
      if (h1 < h0 + minBuildingHeight)
      {
        // warning("Height too small for building " + building.UUID);
        // info("Setting building height to " +
        //     str(minBuildingHeight));
        h1 = h0 + minBuildingHeight;
        numSmallHeights++;
        building.error |= BuildingError::BUILDING_HEIGHT_TOO_LOW;
      }

      // Set building height(s)
      building.Height = h1 - h0;
      building.GroundHeight = h0;

      if (building.Height / sqrt(Geometry::PolygonArea(building.Footprint)) >
          25)
      {
        building.error |= BuildingError::BUILDING_BAD_ASPECT_RATIO;
      }

      // Uncomment for debugging
      // Plotting::Plot(building);
    }

    // Print some statistics
    const size_t n = city.Buildings.size();
    info("Missing ground points for " + str(numMissingGroundPoints) + "/" +
         str(n) + " building(s)");
    info("Missing roof points for " + str(numMissingRoofPoints) + "/" + str(n) +
         " building(s)");
    info("Height too small (adjusted) for " + str(numSmallHeights) + "/" +
         str(n) + " building(s)");
  }

  /// Generate a random city. Used for benchmarking.
  ///
  /// @param city The city
  /// @param dtm Digital Terrain Map
  /// @param numBuildings Number of buildings
  static void
  RandomizeCity(City &city, const GridField &dtm, size_t numBuildings)
  {
    info("Randomizing city...");

    // Some hard-coded dimensions
    const double A = 20.0;  // Maximum building side length
    const double H = 50.0;  // Maximum building height
    const double N = 10000; // Maximum number of attempts

    // Get bounding box of domain
    const BoundingBox2D &bbox = dtm.grid.BoundingBox;
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
          info("Try setting a smaller number of random buildings.");
          error("Unable to randomize city; reached maximum number of "
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
          if (d < 0.1 * A) // Allow big overlaps to get interesting results
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
        city.Buildings.push_back(building);
        centers.push_back(c);

        info("Creating random building " + str(i + 1) + "/" +
             str(numBuildings) + " at c = " + str(c));
        break;
      }
    }
  }
  /// Generate building with given dimensions
  /// @param Poind2D The center of the building footprint
  /// @param a x-width of building
  /// @param b y-width of building
  /// @param height relative (to the ground height) height of the building
  /// @param groundHeight ground height of the building
  /// Note: Absolute building height = height + groundHeight
  static Building GenerateBuilding(
      const Point2D &c, double a, double b, double height, double groundHeight)
  {
    Building building;

    // Set building geometry
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

    // Create ground points and roof points
    const size_t numPoints = 5;
    for (size_t i = 0; i < numPoints; i++)
    {
      building.GroundPoints.push_back(Point3D(c.x, c.y, groundHeight));
      building.RoofPoints.push_back(Point3D(c.x, c.y, groundHeight + height));
    }

    return building;
  }

  static City remove_building_point_outliers_statistical(const City &city,
                                                         size_t neighbours,
                                                         double outlierMargin)
  {
    // Create copy of city
    City _city{city};

    // size_t totalRemoved = 0;
    for (auto &building : _city.Buildings)
    {
      // size_t beforeFilter = building.RoofPoints.size();
      PointCloudProcessor::StatisticalOutlierRemover(building.RoofPoints,
                                                     neighbours, outlierMargin);
      // totalRemoved += (beforeFilter - building.RoofPoints.size());
    }

    return _city;
  }

  static City remove_building_point_outliers_ransac(const City &city,
                                                    double distanceThershold,
                                                    size_t iterations)
  {
    // Create copy of city
    City _city{city};

    // size_t totalRemoved = 0;
    for (auto &building : _city.Buildings)
    {
      // size_t beforeFilter = building.RoofPoints.size();
      PointCloudProcessor::RANSAC_OutlierRemover(building.RoofPoints,
                                                 distanceThershold, iterations);
      // totalRemoved += (beforeFilter - building.RoofPoints.size());
    }

    return _city;
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
};

} // namespace DTCC_BUILDER

#endif
