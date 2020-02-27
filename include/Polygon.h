// Polygon (array of 2D points).
// Copyright (C) 2019 Anders Logg.

#ifndef DTCC_POLYGON_H
#define DTCC_POLYGON_H

#include <cmath>
#include <vector>

#include "Point.h"

namespace DTCC
{

class Polygon
{
public:
  // Array of points
  std::vector<Point2D> Points;

  // Create empty polygon
  Polygon() {}
};

} // namespace DTCC

#endif
