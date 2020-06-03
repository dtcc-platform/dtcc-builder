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
  // Array of vertices
  std::vector<Point2D> Vertices;

  // Create empty polygon
  Polygon() {}
};

} // namespace DTCC

#endif
