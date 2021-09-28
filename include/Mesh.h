// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_H
#define DTCC_MESH_H

#include <vector>

#include "Logging.h"
#include "Point.h"
#include "Simplex.h"
#include "Vector.h"

namespace DTCC
{

/// Mesh2D represents a triangular mesh of a 2D domain.

class Mesh2D : public Printable
{
public:
  /// Array of vertices
  std::vector<Point2D> Vertices{};

  /// Array of cells (triangles)
  std::vector<Simplex2D> Cells{};

  /// Array of cell markers
  std::vector<int> Markers{};

  /// Compute midpoint of cell
  Point2D MidPoint(size_t cellIndex) const
  {
    Vector2D c{};
    c += Vector2D(Vertices[Cells[cellIndex].v0]);
    c += Vector2D(Vertices[Cells[cellIndex].v1]);
    c += Vector2D(Vertices[Cells[cellIndex].v2]);
    c /= 3.0;
    return c;
  }

  /// Pretty-print
  std::string __str__() const override
  {
    return "2D triangular mesh with " + str(Vertices.size()) +
           " vertices and " + str(Cells.size()) + " faces";
  }

  };

  /// Mesh3D represents a triangular mesh of a 3D domain.

  class Mesh3D : public Printable
  {
  public:

    /// Array of vertices
    std::vector<Point3D> Vertices{};

    /// Array of cells (tetrahedra)
    std::vector<Simplex3D> Cells{};

    /// Array of cell markers
    std::vector<int> Markers{};

    /// Compute of cell
    Point3D MidPoint(size_t cellIndex) const
    {
      Vector3D c;
      c += Vector3D(Vertices[Cells[cellIndex].v0]);
      c += Vector3D(Vertices[Cells[cellIndex].v1]);
      c += Vector3D(Vertices[Cells[cellIndex].v2]);
      c += Vector3D(Vertices[Cells[cellIndex].v3]);
      c /= 4.0;
      return c;
    }

    // Pretty-print
    std::string __str__() const override
    {
      return "3D tetrahedral mesh with "
        + str(Vertices.size()) + " vertices and "
        + str(Cells.size()) + " cells";
    }

  };

} // namespace DTCC

#endif
