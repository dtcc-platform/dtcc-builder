// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_SURFACE_H
#define DTCC_SURFACE_H

#include <vector>

#include "Logging.h"
#include "model/Polygon.h"
#include "model/Vector.h"

namespace DTCC_BUILDER
{
class Surface : public Printable
{
public:
  std::vector<Vector3D> vertices{};
  std::vector<std::vector<Vector3D>> holes{};
  Surface() = default;
  virtual ~Surface() {} // make the destructor virtual

  double max_height() const
  {
    double max = -1e9;
    for (auto &v : vertices)
    {
      if (v.z > max)
      {
        max = v.z;
      }
    }
    return max;
  }

  // project to 2D polygon
  Polygon to_polygon() const
  {
    Polygon p;
    for (const auto &v : vertices)
    {
      p.vertices.push_back(Vector2D(v.x, v.y));
    }
    for (const auto &h : holes)
    {
      std::vector<Vector2D> hole;
      for (const auto &v : h)
      {
        hole.push_back(Vector2D(v.x, v.y));
      }
      p.holes.push_back(hole);
    }
    return p;
  }

  /// Pretty-print
  std::string __str__() const override
  {
    return "Surface with " + str(vertices.size()) + " vertices";
  }
};

class MultiSurface
{
public:
  std::vector<Surface> surfaces{};

  MultiSurface() = default;
  virtual ~MultiSurface() {} // make the destructor virtual
};
} // namespace DTCC_BUILDER

#endif