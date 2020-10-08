// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_POLYGON_H
#define DTCC_POLYGON_H

#include <cmath>
#include <vector>

#include "Logging.h"
#include "Point.h"

namespace DTCC
{

class Polygon : public Printable
{
public:
  // Array of vertices
  std::vector<Point2D> Vertices{};

  // Create empty polygon
  Polygon() = default;

  /// Pretty-print
  std::string __str__() const
  {
    std::string s = "[";
    for (size_t i = 0; i < Vertices.size(); i++)
    {
      if (i > 0)
        s += ", ";
      s += str(Vertices[i]);
    }
    s += "]";
    return s;
  }
};

} // namespace DTCC

#endif
