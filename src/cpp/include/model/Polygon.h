// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_POLYGON_H
#define DTCC_POLYGON_H

#include <cmath>
#include <vector>

#include "Logging.h"
#include "Point.h"

namespace DTCC_BUILDER
{

class Polygon : public Printable
{
  typedef std::vector<Point2D> LineString;

public:
  // Array of vertices
  LineString Vertices{};

  std::vector<LineString> Holes{};

  // Create empty polygon
  Polygon() = default;
  virtual ~Polygon() {} // make the destructor virtual

  /// Set new origin (subtract offset)
  void SetOrigin(const Point2D &origin)
  {
    for (auto &p : Vertices)
    {
      p.x -= origin.x;
      p.y -= origin.y;
    }
  }

  /// Pretty-print
  std::string __str__() const override
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

} // namespace DTCC_BUILDER

#endif
