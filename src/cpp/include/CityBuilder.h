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
  /// build city from building footprints, including only building
  /// inside the given bounding box. Note that this does not generate any
  /// building heights, only flat 2D buildings.
  ///
  /// @param city The city
  /// @param footprints Footprints of buildings (polygons)
  /// @param uuids uuids of buildings
  /// @param entity_ids Indices of buildings (in shapefile)
  /// @param bbox Bounding box of domain
  /// @param min_building_distance Minimal distance from building to domain
  /// boundary
  static void build_city(City &city,
                         const std::vector<Polygon> &footprints,
                         const std::vector<std::string> &uuids,
                         const std::vector<int> &entity_ids,
                         const BoundingBox2D &bbox,
                         double min_building_distance,
                         double min_building_size)
  {
    info("Building city...");
    Timer timer("build_city");

    // clear data
    city.buildings.clear();

    // Add buildings
    for (size_t i = 0; i < footprints.size(); i++)
    {
      // Skip if not inside  bounding box
      if (!Geometry::bounding_box_contains_2d(bbox, footprints[i],
                                              min_building_distance))
      {
        warning("Skipping building " + uuids[i] +
                "; outside domain or too close to boundary");
        continue;
      }

      // Add building
      Building building;
      building.footprint = footprints[i];
      building.uuid = uuids[i];
      building.shpfile_id = entity_ids[i];
      if (Geometry::polygon_area(building.footprint) < min_building_size)
      {
        building.error |= BuildingError::BUILDING_TOO_SMALL;
      }

      city.buildings.push_back(building);
    }

    info("Added " + str(city.buildings.size()) + "/" + str(footprints.size()) +
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

    // clear search tree (since it might become invalid)
    city.bbtree.clear();

    // Count some stats
    size_t num_closed = 0;
    size_t num_oriented = 0;
    size_t num_vertex_merged = 0;
    size_t num_edge_merged = 0;
    size_t num_removed = 0;

    // Initialize new city
    City _city;
    _city.name = city.name;

    // Clean buildings
    for (auto &building : city.buildings)
    {
      // Make copy of building
      Building _building{building};

      // Make closed
      num_closed +=
          Polyfix::make_closed(_building.footprint, Constants::epsilon);

      // Make oriented
      num_oriented += Polyfix::make_oriented(_building.footprint);

      // Merge vertices (but skip if only 4 vertices or less)
      if (_building.footprint.vertices.size() > 4)
      {
        num_vertex_merged +=
            Polyfix::merge_vertices(_building.footprint, min_vertex_distance);
      }

      // Merge edges (but skip if only 4 vertices or less)
      if (_building.footprint.vertices.size() > 4)
      {
        num_edge_merged += Polyfix::merge_edges(
            _building.footprint, Constants::footprint_angle_threshold);
      }

      // Keep only valid buildings
      if (_building.valid())
        _city.buildings.push_back(_building);
      else
        num_removed++;
    }

    info("Fixed " + str(num_closed) + "/" + str(city.buildings.size()) +
         " polygons that were not closed");
    info("Fixed " + str(num_oriented) + "/" + str(city.buildings.size()) +
         " polygons that were not oriented");
    info("Merged vertices for " + str(num_vertex_merged) + "/" +
         str(city.buildings.size()) + " polygons");
    info("Merged edges for " + str(num_edge_merged) + "/" +
         str(city.buildings.size()) + " polygons");
    info("Removed " + str(num_removed) + "/" + str(city.buildings.size()) +
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
  /// @param point_cloud Point cloud (unfiltered)
  /// @param ground_margin Margin around building for detecting ground points
  static City compute_building_points(const City &city,
                                      const PointCloud &point_cloud,
                                      double ground_margin,
                                      double ground_outlier_margin)

  {
    info("Computing building points...");
    Timer timer("compute_building_points");

    // Check that point cloud is not empty
    if (point_cloud.points.empty())
      error("empty point cloud");

    // Check that point cloud has classifications
    if (point_cloud.points.size() != point_cloud.classifications.size())
      error("Missing classifications for point cloud");

    bool classified_points =
        (point_cloud.points.size() == point_cloud.classifications.size());
    bool classifed_buildings = false;
    if (classified_points)
      classifed_buildings = point_cloud.has_classification(6);

    auto kdt_timer = Timer("ExtractBuildingPoints: BuildKDTree");
    // build a kd-tree for radius search
    typedef KDTreeVectorOfVectorsAdaptor<std::vector<Vector3D>, double,
                                         2 /* dims */>
        my_kd_tree_t;
    my_kd_tree_t pc_index(2, point_cloud.points, 20 /* max leaf */);
    kdt_timer.stop();

    // Create copy of city
    City _city{city};

    // Iterate over buildings
    for (auto &building : _city.buildings)
    {
      building.ground_points.clear();
      building.roof_points.clear();

      auto centerPoint = Geometry::polygon_center_2d(building.footprint);
      double radius =
          Geometry::polygon_radius_2d(building.footprint, centerPoint);
      radius *= radius;
      radius += ground_margin;

      std::vector<double> query_pt{centerPoint.x, centerPoint.y};
      auto radius_t = Timer("RadiusQuery");
      auto indices_dists = pc_index.radius_query(&query_pt[0], radius);
      radius_t.stop();
      for (auto const &ind_pt : indices_dists)
      {
        size_t idx = ind_pt.first;
        const uint8_t clf = point_cloud.classifications[idx];
        const Vector3D &p_3d = point_cloud.points[idx];
        const Vector2D p_2d{p_3d.x, p_3d.y};

        if (classified_points)
        {
          if (clf == 2 || clf == 9)
          {
            building.ground_points.push_back(p_3d);
          }
          else if (clf == 6 || (!classifed_buildings && clf < 2))
          {
            // auto pc_timer = Timer("PolygoCOntains2D");
            if (Geometry::polygon_contains_2d(building.footprint, p_2d))
            {
              building.roof_points.push_back(p_3d);
            }
            // pc_timer.stop();
          }
        }
        else // unclassified data
        {
          if (Geometry::polygon_contains_2d(building.footprint, p_2d))
          {
            building.roof_points.push_back(p_3d);
          }
          else
          {
            // all points not in the roof polygon are considered ground
            building.ground_points.push_back(p_3d);
          }
        }
      }
    }

    // Remove ground outliers
    size_t num_ground_points = 0;
    size_t num_ground_outliers = 0;
    for (auto &building : _city.buildings)
    {
      // Count total number of points
      num_ground_points += building.ground_points.size();

      // Remove outliers and count total number of outliers
      num_ground_outliers += PointCloudProcessor::remove_outliers(
                                 building.ground_points, ground_outlier_margin)
                                 .size();
    }
    const double outlier_ground_percentage =
        (100.0 * num_ground_outliers) / num_ground_points;
    info("Removed ground point outliers (" + str(outlier_ground_percentage) +
         "%)");

    double pts_pr_sqm;
    double point_coverage;
    size_t too_few = 0;
    for (auto &building : _city.buildings)
    {
      pts_pr_sqm = static_cast<double>(building.roof_points.size()) /
                   Geometry::polygon_area(building.footprint);
      if (pts_pr_sqm < 0.25)
      {
        building.error |= BuildingError::BUILDING_TOO_FEW_POINTS;
        too_few++;
      }
      point_coverage = BuildingProcessor::point_coverage(building, 2.0);
      // info("point_coverage: " + str(point_coverage));
      if (point_coverage < 0.5)
      {
        building.error |= BuildingError::BUILDING_INSUFFICIENT_POINT_COVERAGE;
      }
    }
    info("Number of buildings with too few roof points: " + str(too_few));

    // Sort points by height
    for (auto &building : _city.buildings)
    {
      std::sort(building.ground_points.begin(), building.ground_points.end(),
                [](const Vector3D &p, const Vector3D &q) -> bool
                { return p.z < q.z; });
      std::sort(building.roof_points.begin(), building.roof_points.end(),
                [](const Vector3D &p, const Vector3D &q) -> bool
                { return p.z < q.z; });
    }

    // Compute some statistics
    size_t min_g{std::numeric_limits<size_t>::max()};
    size_t min_r{std::numeric_limits<size_t>::max()};
    size_t max_g{0}, maxR{0};
    size_t sum_g{0}, sumR{0};
    for (const auto &building : _city.buildings)
    {
      // Ground points
      const size_t nG = building.ground_points.size();
      min_g = std::min(min_g, nG);
      max_g = std::max(max_g, nG);
      sum_g += nG;

      // Roof points
      const size_t nR = building.roof_points.size();
      min_r = std::min(min_r, nR);
      maxR = std::max(maxR, nR);
      sumR += nR;
    }
    const double mean_g = static_cast<double>(sum_g) / _city.buildings.size();
    const double mean_r = static_cast<double>(sumR) / _city.buildings.size();

    info("min/mean/max number of ground points per "
         "building is " +
         str(min_g) + "/" + str(mean_g) + "/" + str(max_g));
    info("min/mean/max number of roof points per building "
         "is " +
         str(min_r) + "/" + str(mean_r) + "/" + str(maxR));

    return _city;
  }

  /// Compute heights of buildings from ground and roof points. This
  /// requires that ExtractBuildingPoints() has been called to extract
  /// the points from point cloud data.
  ///
  /// @param city The city
  /// @param dtm Digital Terrain Map, used for ground height if points are
  /// missing
  /// @param ground_percentile Percentile used for setting ground height
  /// @param roof_percentile Percentile used for setting roof height
  static void compute_building_heights(City &city,
                                       const GridField &dtm,
                                       double ground_percentile,
                                       double roof_percentile)
  {
    info("Computing building heights...");
    Timer timer("compute_building_heights");

    // FIXME: Make this a parameter?
    // FIXME: How do we treat this in relation to layer height?
    // FIXME: Also a problem inside MeshGenerator (warning in TrimMesh3D)
    // FIXME: Make this a parameter
    const double min_building_height{2.5};

    // Uncomment for debugging
    // Plotting::Init();

    // Count missing or bad data
    size_t num_missing_ground_points = 0;
    size_t num_missing_roof_points = 0;
    size_t num_small_heights = 0;

    // Iterate over buildings
    for (auto &building : city.buildings)
    {
      // Compute ground height h0
      double h0{0};
      if (building.ground_points.empty())
      {
        // warning("Missing ground points for building " + building.uuid);
        // info("Setting ground height from DTM");
        h0 = dtm(Geometry::polygon_center_2d(building.footprint));
        num_missing_ground_points++;
        building.error |= BuildingError::BUILDING_NO_GROUND_POINTS;
      }
      else
      {
        // Pick percentile from ground points
        sort(building.ground_points.begin(), building.ground_points.end(),
             [](const Vector3D &lhs, const Vector3D &rhs)
             { return lhs.z < rhs.z; });
        h0 = get_percentile(building.ground_points, ground_percentile).z;
      }

      // Compute roof height h1
      double h1{0};
      if (building.roof_points.empty())
      {
        // warning("Missing roof points for building " + building.uuid);
        // info("Setting building height to " +
        //     str(min_building_height) + "m");
        h1 = h0 + min_building_height;
        num_missing_roof_points++;
        building.error |= BuildingError::BUILDING_NO_ROOF_POINTS;
      }
      else
      {
        sort(building.roof_points.begin(), building.roof_points.end(),
             [](const Vector3D &lhs, const Vector3D &rhs)
             { return lhs.z < rhs.z; });
        h1 = get_percentile(building.roof_points, roof_percentile).z;
      }

      // Check that h0 < h1
      if (h1 < h0 + min_building_height)
      {
        // warning("height too small for building " + building.uuid);
        // info("Setting building height to " +
        //     str(min_building_height));
        h1 = h0 + min_building_height;
        num_small_heights++;
        building.error |= BuildingError::BUILDING_HEIGHT_TOO_LOW;
      }

      // Set building height(s)
      building.height = h1 - h0;
      building.ground_height = h0;

      if (building.height / sqrt(Geometry::polygon_area(building.footprint)) >
          25)
      {
        building.error |= BuildingError::BUILDING_BAD_ASPECT_RATIO;
      }

      // Uncomment for debugging
      // Plotting::Plot(building);
    }

    // print some statistics
    const size_t n = city.buildings.size();
    info("Missing ground points for " + str(num_missing_ground_points) + "/" +
         str(n) + " building(s)");
    info("Missing roof points for " + str(num_missing_roof_points) + "/" +
         str(n) + " building(s)");
    info("height too small (adjusted) for " + str(num_small_heights) + "/" +
         str(n) + " building(s)");
  }

  /// Generate a random city. Used for benchmarking.
  ///
  /// @param city The city
  /// @param dtm Digital Terrain Map
  /// @param num_buildings Number of buildings
  static void
  randomize_city(City &city, const GridField &dtm, size_t num_buildings)
  {
    info("Randomizing city...");

    // Some hard-coded dimensions
    const double A = 20.0;  // Maximum building side length
    const double H = 50.0;  // Maximum building height
    const double N = 10000; // Maximum number of attempts

    // Get bounding box of domain
    const BoundingBox2D &bbox = dtm.grid.bounding_box;
    const double dx = bbox.Q.x - bbox.P.x;
    const double dy = bbox.Q.y - bbox.P.y;

    // Iterate over the number of buildings to generate
    std::vector<Vector2D> centers;
    for (size_t i = 0; i < num_buildings; i++)
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
        Vector2D c(bbox.P.x + Utils::random() * dx,
                   bbox.P.y + Utils::random() * dy);

        // Check that we are not too close to other buildings, but
        // note that buildings may actually overlap slightly which may
        // also happen with real-world data.
        bool ok = true;
        for (auto const &p : centers)
        {
          const double d = Geometry::distance_2d(p, c);
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
        const double a = (0.05 + 0.95 * Utils::random()) * A;
        const double b = (0.05 + 0.95 * Utils::random()) * A;
        const double h = (0.25 + 0.75 * Utils::random()) * H;

        // Generate building
        Building building = generate_building(c, a, b, h, dtm(c));

        // Add building
        city.buildings.push_back(building);
        centers.push_back(c);

        info("Creating random building " + str(i + 1) + "/" +
             str(num_buildings) + " at c = " + str(c));
        break;
      }
    }
  }
  /// Generate building with given dimensions
  /// @param Poind2D The center of the building footprint
  /// @param a x-width of building
  /// @param b y-width of building
  /// @param height relative (to the ground height) height of the building
  /// @param ground_height ground height of the building
  /// Note: Absolute building height = height + ground_height
  static Building generate_building(const Vector2D &c,
                                    double a,
                                    double b,
                                    double height,
                                    double ground_height)
  {
    Building building;

    // Set building geometry
    building.footprint.vertices.push_back(
        Vector2D{c.x - 0.5 * a, c.y - 0.5 * b});
    building.footprint.vertices.push_back(
        Vector2D{c.x + 0.5 * a, c.y - 0.5 * b});
    building.footprint.vertices.push_back(
        Vector2D{c.x + 0.5 * a, c.y + 0.5 * b});
    building.footprint.vertices.push_back(
        Vector2D{c.x - 0.5 * a, c.y + 0.5 * b});
    building.height = height;
    building.ground_height = ground_height;

    // Create ground points and roof points
    const size_t num_points = 5;
    for (size_t i = 0; i < num_points; i++)
    {
      building.ground_points.push_back(Vector3D(c.x, c.y, ground_height));
      building.roof_points.push_back(
          Vector3D(c.x, c.y, ground_height + height));
    }

    return building;
  }

  static City remove_building_point_outliers_statistical(const City &city,
                                                         size_t neighbours,
                                                         double outlier_margin)
  {
    // Create copy of city
    City _city{city};

    // size_t totalRemoved = 0;
    for (auto &building : _city.buildings)
    {
      // size_t beforeFilter = building.roof_points.size();
      PointCloudProcessor::statistical_outlier_remover(
          building.roof_points, neighbours, outlier_margin);
      // totalRemoved += (beforeFilter - building.roof_points.size());
    }

    return _city;
  }

  static City remove_building_point_outliers_ransac(const City &city,
                                                    double distance_thershold,
                                                    size_t iterations)
  {
    // Create copy of city
    City _city{city};

    // size_t totalRemoved = 0;
    for (auto &building : _city.buildings)
    {
      // size_t beforeFilter = building.roof_points.size();
      PointCloudProcessor::ransac_outlier_remover(
          building.roof_points, distance_thershold, iterations);
      // totalRemoved += (beforeFilter - building.roof_points.size());
    }

    return _city;
  }

private:
  // Get percentile object from array. It is assumed that the array is ordered.
  template <class T>
  static T get_percentile(const std::vector<T> &array, double percentile)
  {
    size_t index = std::max(0.0, percentile * array.size());
    index = std::min(index, array.size() - 1);
    return array[index];
  }
};

} // namespace DTCC_BUILDER

#endif
