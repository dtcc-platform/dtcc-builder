// Polygon (array of 2D points).
// Copyright (C) 2019 Anders Logg.

#ifndef VC_POLYGON_H
#define VC_POLYGON_H

#include <cmath>
#include <vector>

#include "Point.h"

namespace VirtualCity
{

class Polygon
{
public:
  // Array of points
  std::vector<Point2D> Points;

  // Create empty polygon
  Polygon() {}
};

} // namespace VirtualCity

#endif
