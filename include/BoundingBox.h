// Axis aligned bounding box
// Copyright (C) 2020 Anders Logg

#ifndef DTCC_BOUNDING_BOX_H
#define DTCC_BOUNDING_BOX_H

#include <limits>
#include "Point.h"

namespace DTCC
{

  class BoundingBox2D
  {
  public:
    // First ("lower left") corner
    Point2D P;

    // Second ("upper right") corner
    Point2D Q;

    // Create empty bounding box
    BoundingBox2D() {}

    // Create bounding box of polygon
    BoundingBox2D(const Polygon& polygon)
    {
      constexpr double max = std::numeric_limits<double>::max();
      P.x = P.y = max;
      Q.x = Q.y = -max;
      for (const auto& p: polygon.Points)
      {
        P.x = std::min(P.x, p.x);
        P.y = std::min(P.y, p.y);
        Q.x = std::max(Q.x, p.x);
        Q.y = std::max(Q.y, p.y);
      }
    }
  };

  class BoundingBox3D
  {
  public:
    // First ("lower left") corner
    Point3D P;

    // Second ("upper right") corner
    Point3D Q;

    // Create empty bounding box
    BoundingBox3D() {}
  };

}

#endif
