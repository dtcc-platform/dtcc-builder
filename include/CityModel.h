// Representation of a 2.5D city model.
// Copyright (C) 2019 Anders Logg.

#ifndef DTCC_CITY_MODEL_H
#define DTCC_CITY_MODEL_H

#include <iomanip>
#include <string>
#include <vector>

#include "Building.h"
#include "Geometry.h"
#include "Point.h"
#include "BoundingBoxTree.h"

namespace DTCC
{

class CityModel
{
public:
  // List of buildings
  std::vector<Building> Buildings;

  // Build search tree
  void BuildSearchTree() const
  {
    // Create bounding boxes for all building footprints
    std::vector<BoundingBox2D> bboxes;
    for (const auto& building: Buildings)
    {
      BoundingBox2D bbox(building.Footprint);
      bboxes.push_back(bbox);
    }

    // Build bounding box tree
    bbtree.Build(bboxes);
    std::cout << bbtree << std::endl;
  }

  // Find building containing point (inside footprint), returning -1
  // if the point is not inside any building.
  int FindBuilding(const Point2D &p) const
  {
    // Check that search tree has been created
    if (bbtree.Nodes.size() == 0)
      throw std::runtime_error("Missing search tree; call BuildSearchTree()");

    // Find candidate buildings from search tree
    std::vector<size_t> indices = bbtree.Find(p);

    // Check candidate buildings
    for (const auto index: indices)
    {
      if (Geometry::PolygonContains2D(Buildings[index].Footprint, p))
        return index;
    }

    // Point not inside a building
    return -1;
  }

  // Compute center of city model
  Point2D Center() const
  {
    Point2D c;
    size_t numPoints = 0;
    for (auto const &building : Buildings)
    {
      for (auto const &p : building.Footprint.Points)
      {
        c += p;
        numPoints += 1;
      }
    }
    c /= numPoints;
    return c;
  }

  // Compute radius of city model (relative to center)
  double Radius(const Point2D &center) const
  {
    double r2max = 0.0;
    for (auto const &building : Buildings)
    {
      for (auto const &p : building.Footprint.Points)
      {
        const double r2 = Geometry::SquaredDistance2D(p, center);
        if (r2 > r2max)
          r2max = r2;
      }
    }
    return std::sqrt(r2max);
  }

private:

  // Bounding box tree (used by find building)
  mutable BoundingBoxTree2D bbtree;

};

std::ostream &operator<<(std::ostream &stream, const CityModel &cityModel)
{
  const Point2D c = cityModel.Center();
  const double r = cityModel.Radius(c);
  stream << "CityModel with " << cityModel.Buildings.size()
         << " buildings and radius R = " << std::fixed << std::setprecision(2)
         << r;
  return stream;
}

} // namespace DTCC

#endif
