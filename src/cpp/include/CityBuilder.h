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
#include "GEOS.h"
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
    info("CityBuilder: Building city...");
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

    info("CityBuilder: Added " + str(city.Buildings.size()) + "/" +
         str(footprints.size()) + " buildings inside bounding box");
  }

  /// Clean city by making sure that all building footprints
  /// are closed and counter-clockwise oriented.
  ///
  /// @param city The city
  /// @param minimalVertexDistance Minimal vertex distance
  ///
  /// Developer note: This may be optimized by avoiding copying
  static void CleanCity(City &city, double minimalVertexDistance)
  {
    info("CityBuilder: Cleaning city...");
    Timer timer("CleanCity");

    // Clear search tree (since it might become invalid)
    city.bbtree.Clear();

    // Count some stats
    size_t numClosed = 0;
    size_t numOriented = 0;
    size_t numVertexMerged = 0;
    size_t numEdgeMerged = 0;
    size_t numRemoved = 0;

    // Clean buildings
    for (auto &building : city.Buildings)
    {
      // Make closed
      numClosed += Polyfix::MakeClosed(building.Footprint, Constants::Epsilon);

      // Make oriented
      numOriented += Polyfix::MakeOriented(building.Footprint);

      // Merge vertices (but skip if only 4 vertices or less)
      if (building.Footprint.Vertices.size() > 4)
      {
        numVertexMerged +=
            Polyfix::MergeVertices(building.Footprint, minimalVertexDistance);
      }

      // Merge edges (but skip if only 4 vertices or less)
      if (building.Footprint.Vertices.size() > 4)
      {
        numEdgeMerged += Polyfix::MergeEdges(
            building.Footprint, Constants::FootprintAngleThreshold);
      }
    }

    // Keep only valid buildings
    std::vector<Building> _buildings{city.Buildings};
    city.Buildings.clear();
    for (auto &building : _buildings)
    {
      if (building.Valid())
        city.Buildings.push_back(building);
      else
        numRemoved++;
    }

    info("CityBuilder: Fixed " + str(numClosed) + "/" +
         str(city.Buildings.size()) + " polygons that were not closed");
    info("CityBuilder: Fixed " + str(numOriented) + "/" +
         str(city.Buildings.size()) + " polygons that were not oriented");
    info("CityBuilder: Merged vertices for " + str(numVertexMerged) + "/" +
         str(city.Buildings.size()) + " polygons");
    info("CityBuilder: Merged edges for " + str(numEdgeMerged) + "/" +
         str(city.Buildings.size()) + " polygons");
    info("CityBuilder: Removed " + str(numRemoved) + "/" +
         str(city.Buildings.size()) +
         " buildings (invalid/too small after cleaning)");
  }

  /// Simplify city by merging all buildings that are closer than
  /// a given distance. When merging buildings, the number of buildings
  /// will decrease. Ground points and roof points are also merged and
  /// heights are set to min/max values of the merged buildings.
  static void SimplifyCity(City &city,
                           const BoundingBox2D &bbox,
                           double minimalBuildingDistance,
                           double minimalVertexDistance)
  {
    info("CityBuilder: Simplifying city...");
    Timer timer("SimplifyCity");

    // Clear search tree (since it might become invalid)
    city.bbtree.Clear();

    // Merge buildings if too close
    MergeCity(city, bbox, minimalBuildingDistance);
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
  /// @param city The city
  /// @param pointCloud Point cloud (unfiltered)
  /// @param groundMargin Margin around building for detecting ground points
  static void ExtractBuildingPoints(City &city,
                                    const PointCloud &pointCloud,
                                    double groundMargin,
                                    double groundOutlierMargin)

  {
    info("CityBuilder: Extracting building points...");
    Timer timer("ExtractBuildingPoints");

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
    for (auto &building : city.Buildings)
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
    for (auto &building : city.Buildings)
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
    info("CityBuilder: Removed ground point outliers (" +
         str(outlierGroundPercentage) + "%)");

    double ptsPrSqm;
    double pointCoverage;
    size_t tooFew = 0;
    for (auto &building : city.Buildings)
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
    info("CityBuilder: Number of buildings with too few roof points: " +
         str(tooFew));

    // Sort points by height
    for (auto &building : city.Buildings)
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
    for (const auto &building : city.Buildings)
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
    const double meanG = static_cast<double>(sumG) / city.Buildings.size();
    const double meanR = static_cast<double>(sumR) / city.Buildings.size();

    info("CityBuilder: min/mean/max number of ground points per "
         "building is " +
         str(minG) + "/" + str(meanG) + "/" + str(maxG));
    info("CityBuilder: min/mean/max number of roof points per building "
         "is " +
         str(minR) + "/" + str(meanR) + "/" + str(maxR));
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
    info("CityBuilder: Computing building heights...");
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
        // info("CityBuilder: Setting ground height from DTM");
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
        // info("CityBuilder: Setting building height to " +
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
        // info("CityBuilder: Setting building height to " +
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
    info("CityBuilder: Missing ground points for " +
         str(numMissingGroundPoints) + "/" + str(n) + " building(s)");
    info("CityBuilder: Missing roof points for " + str(numMissingRoofPoints) +
         "/" + str(n) + " building(s)");
    info("CityBuilder: Height too small (adjusted) for " +
         str(numSmallHeights) + "/" + str(n) + " building(s)");
  }

  /// Generate a random city. Used for benchmarking.
  ///
  /// @param city The city
  /// @param dtm Digital Terrain Map
  /// @param numBuildings Number of buildings
  static void
  RandomizeCity(City &city, const GridField &dtm, size_t numBuildings)
  {
    info("CityBuilder: Randomizing city...");

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

  static void BuildingPointsOutlierRemover(City &city,
                                           size_t neighbours,
                                           double outlierMargin,
                                           bool verbose = false)
  {
    if (verbose)
      info("CityBuilder: BuildingPointsOutlierRemover");
    Timer("BuildingPointsOutlierRemover");
    size_t totalRemoved = 0;
    for (auto &building : city.Buildings)
    {
      size_t beforeFilter = building.RoofPoints.size();
      PointCloudProcessor::StatisticalOutlierRemover(
          building.RoofPoints, neighbours, outlierMargin, verbose);
      totalRemoved += (beforeFilter - building.RoofPoints.size());
    }
    if (verbose)
    {
      info("BuildingPointsOutlierRemove filtered a total of " +
           str(totalRemoved) + " points from  " + str(city.Buildings.size()) +
           " buildings");
    }
  }

  static void BuildingPointsRANSACOutlierRemover(City &city,
                                                 double distanceThershold,
                                                 size_t iterations,
                                                 bool verbose = false)
  {
    if (verbose)
      info("CityBuilder: BuildingPointsRANSACOutlierRemover");
    Timer("BuildingPointsRANSACOutlierRemover");
    size_t totalRemoved = 0;
    for (auto &building : city.Buildings)
    {
      size_t beforeFilter = building.RoofPoints.size();
      PointCloudProcessor::RANSAC_OutlierRemover(building.RoofPoints,
                                                 distanceThershold, iterations);
      totalRemoved += (beforeFilter - building.RoofPoints.size());
    }
    info("BuildingPointsRANSACOutlierRemover remove " + str(totalRemoved) +
         " points");
  }

private:
  // Merge all buildings closer than a given distance
  static void MergeCity(City &city,
                        const BoundingBox2D &bbox,
                        double minimalBuildingDistance)
  {
    info("CityBuilder: Merging buildings...");

    // Initialize GEOS
    GEOS::Init();

    // Avoid using sqrt for efficiency
    const double tol2 = minimalBuildingDistance * minimalBuildingDistance;

    // Get buildings
    std::vector<Building> &buildings = city.Buildings;

    // Counters
    size_t numMerged = 0;
    size_t numCompared = 0;

    // Initialize grid
    // Note: Grid size needs to be *at least* minimal distance
    // Note: Factor 4 seems to be a good choice (tested using dtcc-bench-run)
    double h = 4.0 * ComputeMeanBuildingSize(buildings);
    h = std::max(h, minimalBuildingDistance + Constants::Epsilon);
    size_t nX = static_cast<size_t>((bbox.Q.x - bbox.P.x) / h) + 1;
    size_t nY = static_cast<size_t>((bbox.Q.y - bbox.P.y) / h) + 1;

    if (nX <= 1 || nY <= 1)
    {
      // needed for meshing small areas
      nX = static_cast<size_t>((bbox.Q.x - bbox.P.x) / (h / 4)) + 1;
      nY = static_cast<size_t>((bbox.Q.y - bbox.P.y) / (h / 4)) + 1;
    }

    Grid grid(bbox, nX, nY);

    // Initialize bins
    std::vector<std::unordered_set<size_t>> building2bins{buildings.size()};
    std::vector<std::unordered_set<size_t>> bin2buildings{grid.NumVertices()};
    for (size_t i = 0; i < buildings.size(); i++)
      UpdateBinning(building2bins, bin2buildings, i, buildings[i], grid);

    // Create queue of indices to check
    std::queue<size_t> indices{};
    for (size_t i = 0; i < buildings.size(); i++)
      indices.push(i);

    // Process queue until empty
    while (!indices.empty())
    {
      // Pop index of next building to check
      const size_t i = indices.front();
      indices.pop();

      // Get neighbor indices
      std::unordered_set<size_t> neighbors{
          GetNeighbors(i, building2bins, bin2buildings)};

      // Iterate over neighbors
      for (size_t j : neighbors)
      {
        // Skip building itself
        if (i == j)
          continue;

        // Skip if merged with other building (size set to 0)
        if (buildings[j].Empty())
          continue;

        // Compute distance
        const Polygon &Pi = buildings[i].Footprint;
        const Polygon &Pj = buildings[j].Footprint;
        const double d2 = Geometry::SquaredDistance2D(Pi, Pj);
        numCompared++;

        // Merge if distance is small
        if (d2 < tol2)
        {
          debug("CityBuilder: Buildings " + str(i) + " and " + str(j) +
                " are too close, merging");

          // Merge buildings
          buildings[i].AttachedUUIDs.push_back(buildings[j].UUID);
          MergeBuildings(buildings[i], buildings[j], minimalBuildingDistance);
          numMerged++;

          // Update binning
          UpdateBinning(building2bins, bin2buildings, i, buildings[i], grid);

          // Add building back to queue
          indices.push(i);
        }
      }
    }

    // Extract non-empty polygons (might be done more efficiently)
    std::vector<Building> mergedBuildings;
    for (const auto &building : buildings)
    {
      if (building.Valid())
        mergedBuildings.push_back(building);
      else if (!building.Empty())
        warning("Building " + building.UUID +
                " has non-empty footprint but less than 3 vertices, skipping");
    }

    // Overwrite buildings
    city.Buildings = mergedBuildings;

    // Finish GEOS
    GEOS::Finish();

    info("CityBuilder: " + str(numMerged) + " building pair(s) were merged");
    info("CityBuilder: " + str(numCompared) +
         " pair(s) of buildings were checked");
  }

  // Compute mean building size (from bounding boxes)
  static double ComputeMeanBuildingSize(std::vector<Building> &buildings)
  {
    double meanBuildingSize = 0.0;
    for (const auto &building : buildings)
    {
      BoundingBox2D bbox(building.Footprint.Vertices);
      meanBuildingSize += std::max(bbox.Q.x - bbox.P.x, bbox.Q.y - bbox.P.y);
    }
    meanBuildingSize /= static_cast<double>(buildings.size());
    return meanBuildingSize;
  }

  // Update binning for for building
  static void
  UpdateBinning(std::vector<std::unordered_set<size_t>> &building2bins,
                std::vector<std::unordered_set<size_t>> &bin2buildings,
                size_t buildingIndex,
                const Building &building,
                const Grid &grid)
  {
    // Compute bounding box of building
    BoundingBox2D bbox(building.Footprint.Vertices);

    // Get grid cell size
    const double hx = grid.XStep;
    const double hy = grid.YStep;

    // Get grid indices for bounding box
    long int ixMin{}, iyMin{};
    long int ixMax{}, iyMax{};
    grid.Point2Index(ixMin, iyMin, bbox.P);
    grid.Point2Index(ixMax, iyMax, bbox.Q);

    // Check margin
    double xMin = grid.BoundingBox.P.x + ixMin * hx;
    double yMin = grid.BoundingBox.P.y + iyMin * hy;
    double xMax = grid.BoundingBox.P.x + ixMax * hx;
    double yMax = grid.BoundingBox.P.y + iyMax * hy;
    if (xMin - bbox.P.x + Constants::Epsilon > 0.0)
      ixMin -= 1;
    if (yMin - bbox.P.y + Constants::Epsilon > 0.0)
      iyMin -= 1;
    if (bbox.Q.x - xMax + Constants::Epsilon > 0.0)
      ixMax += 1;
    if (bbox.Q.y - yMax + Constants::Epsilon > 0.0)
      iyMax += 1;

    // Check overflow
    if (ixMin < 0)
      ixMin = 0;
    if (iyMin < 0)
      iyMin = 0;
    if (ixMax >= grid.XSize)
      ixMax = grid.XSize - 1;
    if (iyMax >= grid.YSize)
      iyMax = grid.YSize - 1;

    // Add to bins
    for (long int ix = ixMin; ix <= ixMax; ix++)
    {
      for (long int iy = iyMin; iy <= iyMax; iy++)
      {
        const long int binIndex = grid.Index2Index(ix, iy);
        building2bins[buildingIndex].insert(binIndex);
        bin2buildings[binIndex].insert(buildingIndex);
      }
    }

    // Sanity check: These numbers should never be larger
    // than 0 and only rarely smaller than -0.5
    const long int minIndex = grid.Index2Index(ixMin, iyMin);
    const long int maxIndex = grid.Index2Index(ixMax, iyMax);
    const Point2D P = grid.Index2Point(minIndex);
    const Point2D Q = grid.Index2Point(maxIndex);
    const double dxMin = (P.x - bbox.P.x) / hx;
    const double dxMax = (bbox.Q.x - Q.x) / hx;
    const double dyMin = (P.y - bbox.P.y) / hy;
    const double dyMax = (bbox.Q.y - Q.y) / hy;
    // std::cout << "CHECK: " << dxMin << " " << dxMax << " " << dyMin << " " <<
    //  dyMax << std::endl;
    assert(dxMin < 0.0);
    assert(dxMax < 0.0);
    assert(dyMin < 0.0);
    assert(dyMax < 0.0);
  }

  // Get neighbors of building (buildings with overlapping bins)
  static std::unordered_set<size_t>
  GetNeighbors(size_t buildingIndex,
               const std::vector<std::unordered_set<size_t>> &building2bins,
               const std::vector<std::unordered_set<size_t>> &bin2buildings)
  {
    std::unordered_set<size_t> indices{};
    for (const auto binIndex : building2bins[buildingIndex])
      for (const auto index : bin2buildings[binIndex])
        indices.insert(index);
    return indices;
  }

  // Merge two buildings, replacing the first building and clearing the second.
  // Note that we don't yet handle UUID and ID (using value for first building).
  static void MergeBuildings(Building &building0,
                             Building &building1,
                             double minimalBuildingDistance)
  {
    // Compute merged polygon (old implementation)
    // Polygon mergedPolygon = Polyfix::MergePolygons(
    //  building0.Footprint, building1.Footprint, minimalBuildingDistance);

    // Compute merged polygon (using GEOS)
    Polygon mergedPolygon = GEOS::MergePolygons(
        building0.Footprint, building1.Footprint, minimalBuildingDistance);

    // Set merged polygon
    building0.Footprint = mergedPolygon;

    // Set merged ground and roof points (just append)
    for (const auto &p : building1.GroundPoints)
      building0.GroundPoints.push_back(p);
    for (const auto &p : building1.RoofPoints)
      building0.RoofPoints.push_back(p);

    // Set merged heights (use min and max)
    const double h0 = std::min(building0.MinHeight(), building1.MinHeight());
    const double h1 = std::max(building0.MaxHeight(), building1.MaxHeight());
    building0.GroundHeight = h0;
    building0.Height = h1 - h0;

    // Erase second building
    building1.Clear();
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
        if (d < Constants::Epsilon)
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
