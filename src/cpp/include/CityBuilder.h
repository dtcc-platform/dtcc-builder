// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_CITY_BUILDER_H
#define DTCC_CITY_BUILDER_H

#include <algorithm>
#include <iostream>
#include <map>
#include <mutex>
#include <queue>
#include <set>
#include <thread>
#include <unordered_set>
#include <vector>

// #ifdef _OPENMP
// #include <omp.h>
// #endif

#include "BuildingProcessor.h"
#include "CityProcessor.h"
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

  /// Extract roof points from point cloud.
  ///
  ///
  /// The roof points of a building are defined as all points
  /// of class 6 (Building) that fall within the building
  /// footprint. However, since that classification seems to
  /// be missing in the data from LM, we are currently using
  /// all points (except class 2 and 9).
  ///
  /// @param city The city
  /// @param point_cloud Point cloud (unfiltered)
  static std::vector<std::vector<Vector3D>>
  extract_building_points(const std::vector<Polygon> &footprints,
                          const std::vector<Vector3D> &points)
  {
    info("Computing building points...");
    Timer timer("compute_building_points");

    // Check that point cloud is not empty
    if (points.empty())
      error("empty point cloud");

    // auto tile_city_timer = Timer("Tile City");
    // auto tiles_city = CityProcessor::tile_citymodel(
    //     city, point_cloud, point_cloud.bounding_box, 4, 4);
    // tile_city_timer.stop();
    auto kdt_timer = Timer("ExtractBuildingPoints: BuildKDTree");
    // build a kd-tree for radius search
    typedef KDTreeVectorOfVectorsAdaptor<std::vector<Vector3D>, double,
                                         2 /* dims */>
        my_kd_tree_t;
    my_kd_tree_t pc_index(2, points, 20 /* max leaf */);
    kdt_timer.stop();

    // Iterate over buildings
    std::vector<std::vector<Vector3D>> building_points;
    for (auto &footprint : footprints)
    {
      std::vector<Vector3D> roof_points;
      auto centerPoint = Geometry::polygon_center_2d(footprint);
      double radius = Geometry::polygon_radius_2d(footprint, centerPoint);
      radius *= radius;

      std::vector<double> query_pt{centerPoint.x, centerPoint.y};
      auto radius_t = Timer("RadiusQuery");
      auto indices_dists = pc_index.radius_query(&query_pt[0], radius);
      radius_t.stop();
      for (auto const &ind_pt : indices_dists)
      {
        size_t idx = ind_pt.first;
        const Vector3D &p_3d = points[idx];
        const Vector2D p_2d{p_3d.x, p_3d.y};
        if (Geometry::polygon_contains_2d(footprint, p_2d))
        {
          roof_points.push_back(p_3d);
        }
      }
      building_points.push_back(roof_points);
    }
    return building_points;
  }

  // static City compute_building_points_parallel(const City &city,
  //                                              const PointCloud &point_cloud,
  //                                              size_t x_tiles,
  //                                              size_t y_tiles)
  // {
  //   info("Computing building points...");
  //   Timer timer("compute_building_points in parallel");
  //   Timer tile_timer("Tile City");
  //   auto tiles = CityProcessor::tile_citymodel(
  //       city, point_cloud, point_cloud.bounding_box, x_tiles, y_tiles);
  //   tile_timer.stop();
  //   City out_city;
  //   std::mutex out_city_mutex;
  //   std::vector<std::thread> threads;
  //   for (auto &tile : tiles)
  //   {
  //     if (tile.first.buildings.empty())
  //       continue;
  //     threads.emplace_back(
  //         [&]
  //         {
  //           auto tile_city = compute_building_points(tile.first,
  //           tile.second); std::lock_guard<std::mutex> lock(out_city_mutex);
  //           out_city.merge(tile_city);
  //         });
  //   }
  //   std::for_each(threads.begin(), threads.end(),
  //                 [](std::thread &t) { t.join(); });
  //   timer.stop();
  //   // Timer::report("compute_building_points in parallel");
  //   return out_city;
  // }

  //

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

    // #pragma omp parallel for
    for (size_t i = 0; i < _city.buildings.size(); i++)
    {
      auto &building = _city.buildings[i];
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
    // #pragma omp parallel for
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
