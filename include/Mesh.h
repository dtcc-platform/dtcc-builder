// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_H
#define DTCC_MESH_H

#include "Point.h"
#include "Simplex.h"
#include <vector>

namespace DTCC
{

class Mesh2D
{
public:
  /// Array of points (vertices)
  std::vector<Point2D> Points;

  /// Array of cells (triangles)
  std::vector<Simplex2D> Cells;

  /// Array of domain markers
  std::vector<int> DomainMarkers;

  /// Compute cell midpoint
  Point2D MidPoint(const Simplex2D &Cell) const
  {
    Point2D c;
    c += Points[Cell.v0];
    c += Points[Cell.v1];
    c += Points[Cell.v2];
    c /= 3.0;
    return c;
  }
};

class Mesh3D
{
public:
  /// Array of points (vertices)
  std::vector<Point3D> Points;

  /// Array of cells (tetrahedra)
  std::vector<Simplex3D> Cells;

  /// Array of domain markers
  std::vector<int> DomainMarkers;

  /// Compute cell midpoint
  Point3D MidPoint(const Simplex3D &Cell) const
  {
    Point3D c;
    c += Points[Cell.v0];
    c += Points[Cell.v1];
    c += Points[Cell.v2];
    c += Points[Cell.v3];
    c /= 4.0;
    return c;
  }
};

std::ostream &operator<<(std::ostream &stream, const Mesh2D &m)
{
  stream << "2D mesh with " << m.Points.size() << " points and "
         << m.Cells.size() << " cells (triangles)";
  return stream;
}

std::ostream &operator<<(std::ostream &stream, const Mesh3D &m)
{
  stream << "3D mesh with " << m.Points.size() << " points and "
         << m.Cells.size() << " cells (tetrahedra)";
  return stream;
}

} // namespace DTCC

#endif
