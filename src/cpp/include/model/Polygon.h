// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_POLYGON_H
#define DTCC_POLYGON_H

#include <cmath>
#include <vector>

#include "Logging.h"
#include "Vector.h"

namespace DTCC_BUILDER
{

class Polygon : public Printable
{
  typedef std::vector<Vector2D> LineString;

public:
  // Array of vertices
  LineString vertices{};

  std::vector<LineString> holes{};

  // Create empty polygon
  Polygon() = default;
  virtual ~Polygon() {} // make the destructor virtual

  /// Set new origin (subtract offset)
  void set_origin(const Vector2D &origin)
  {
    for (auto &p : vertices)
    {
      p.x -= origin.x;
      p.y -= origin.y;
    }
  }

  /// Pretty-print
  std::string __str__() const override
  {
    std::string s = "[";
    for (size_t i = 0; i < vertices.size(); i++)
    {
      if (i > 0)
        s += ", ";
      s += str(vertices[i]);
    }
    s += "]";
    return s;
  }
};

} // namespace DTCC_BUILDER

#endif
