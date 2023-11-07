// Copyright (C) 2021 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_CITY_MODEL_PROCESSOR_H
#define DTCC_CITY_MODEL_PROCESSOR_H

#include <vector>

#include "Geometry.h"

#include "BoundingBox.h"
#include "BoundingBoxTree.h"
#include "model/Building.h"
#include "model/City.h"
#include "model/PointCloud.h"
namespace DTCC_BUILDER
{

class CityProcessor
{
public:
  static std::vector<std::pair<City, PointCloud>>
  tile_citymodel(const City &city,
                 const PointCloud &pc,
                 const BoundingBox2D &bounds,
                 size_t x_tiles,
                 size_t y_tiles)
  {
    std::vector<std::pair<City, PointCloud>> tiled_citymodels;
    auto tiles = tile_bounding_box(bounds, x_tiles, y_tiles);
    for (size_t idx = 0; idx < tiles.size(); idx++)
    {
      auto city = City();
      city.name = "tile_" + std::to_string(idx);
      tiled_citymodels.push_back(std::make_pair(city, PointCloud()));
    }
    BoundingBoxTree2D tile_tree;
    tile_tree.build(tiles);
    city.build_search_tree();

    auto building_tile_indices = city.bbtree.find(tile_tree);
    std::unordered_map<size_t, size_t> building_index_to_tile_index;
    for (size_t i = 0; i < building_tile_indices.size(); i++)
    {
      building_index_to_tile_index[building_tile_indices[i].first] =
          building_tile_indices[i].second;
    }
    // info("Found " + str(building_tile_indices.size()) + " buildings in
    // tiles");

    for (size_t i = 0; i < city.buildings.size(); i++)
    {
      auto building = city.buildings[i];
      auto tile_index = building_index_to_tile_index[i];
      tiled_citymodels[tile_index].first.buildings.push_back(building);
    }

    std::vector<BoundingBox2D> tile_buildings_bounds;
    for (size_t i = 0; i < tiled_citymodels.size(); i++)
    {
      auto buildings = tiled_citymodels[i].first.buildings;
      auto tile_bb = BoundingBox2D(buildings[0].footprint);
      for (size_t j = 1; j < buildings.size(); j++)
      {
        tile_bb.union_with(BoundingBox2D(buildings[j].footprint));
      }
      tile_buildings_bounds.push_back(tile_bb);
    }
    BoundingBoxTree2D tile_buildings_bounds_tree;
    tile_buildings_bounds_tree.build(tile_buildings_bounds);
    for (size_t idx = 0; idx < pc.points.size(); idx++)
    {
      auto point = pc.points[idx];
      auto cls = pc.classifications[idx];
      auto tile_index = tile_buildings_bounds_tree.find(point);
      for (auto const &i : tile_index)
      {
        tiled_citymodels[i].second.points.push_back(point);
        tiled_citymodels[i].second.classifications.push_back(cls);
      }
    }

    return tiled_citymodels;
  }

private:
  static std::vector<BoundingBox2D>
  tile_bounding_box(const BoundingBox2D &bounds, size_t x_tiles, size_t y_tiles)
  {
    std::vector<BoundingBox2D> tiles;
    double x_step = (bounds.Q.x - bounds.P.x) / x_tiles;
    double y_step = (bounds.Q.y - bounds.P.y) / y_tiles;
    for (size_t x = 0; x < x_tiles; x++)
    {
      for (size_t y = 0; y < y_tiles; y++)
      {
        BoundingBox2D tile_bb;
        tile_bb.P.x = bounds.P.x + x * x_step;
        tile_bb.Q.x = bounds.P.x + (x + 1) * x_step;

        tile_bb.P.y = bounds.P.y + y * y_step;
        tile_bb.Q.y = bounds.P.y + (y + 1) * y_step;

        tiles.push_back(tile_bb);
      }
    }
    return tiles;
  }
};

} // namespace DTCC_BUILDER
#endif
