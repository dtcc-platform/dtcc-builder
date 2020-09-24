// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_H
#define DTCC_MESH_H

#include <vector>

#include "Vector.h"
#include "Simplex.h"
#include "Logging.h"

namespace DTCC
{

  class Mesh2D : public Printable
  {
  public:
    /// Array of vertices
    std::vector<Vector2D> Vertices;

    /// Array of cells (triangles)
    std::vector<Simplex2D> Cells;

    /// Array of domain markers
    std::vector<int> DomainMarkers;

    /// Compute cell midpoint
    Vector2D MidPoint(const Simplex2D &cell) const
    {
      Vector2D c;
      c += Vertices[cell.v0];
      c += Vertices[cell.v1];
      c += Vertices[cell.v2];
      c /= 3.0;
      return c;
    }

    /// Pretty-print
    std::string __str__() const
    {
      return "2D triangular mesh with "
        + str(Vertices.size()) + " vertices and "
        + str(Cells.size()) + " faces";
    }

  };

  class Mesh3D : public Printable
  {
  public:

    /// Array of vertices
    std::vector<Vector3D> Vertices;

    /// Array of cells (tetrahedra)
    std::vector<Simplex3D> Cells;

    /// Array of domain markers
    std::vector<int> DomainMarkers;

    /// Compute cell midpoint
    Vector3D MidPoint(const Simplex3D &Cell) const
    {
      Vector3D c;
      c += Vertices[Cell.v0];
      c += Vertices[Cell.v1];
      c += Vertices[Cell.v2];
      c += Vertices[Cell.v3];
      c /= 4.0;
      return c;
    }

    // Pretty-print
    std::string __str__() const
    {
      return "3D tetrahedral mesh with "
        + str(Vertices.size()) + " vertices and "
        + str(Cells.size()) + " cells";
    }

  };

} // namespace DTCC

#endif
