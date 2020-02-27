// Surface classes for 2D and 3D (boundaries of 2D and 3D meshes).
// Copyright (C) 2019 Anders Logg.

#ifndef DTCC_SURFACE_H
#define DTCC_SURFACE_H

#include "Point.h"
#include "Simplex.h"
#include <vector>

namespace DTCC
{

class Surface2D
{
public:
  // List of points (vertices)
  std::vector<Point2D> Points;

  // List of cells (segments)
  std::vector<Simplex1D> Cells;
};

class Surface3D
{
public:
  // List of points (vertices)
  std::vector<Point3D> Points;

  // List of cells (triangles)
  std::vector<Simplex2D> Cells;
};

std::ostream &operator<<(std::ostream &stream, const Surface2D &s)
{
  stream << "2D surface with " << s.Points.size() << " points and "
         << s.Cells.size() << " cells (segments)";
  return stream;
}

std::ostream &operator<<(std::ostream &stream, const Surface3D &s)
{
  stream << "3D surface with " << s.Points.size() << " points and "
         << s.Cells.size() << " cells (triangles)";
  return stream;
}

} // namespace DTCC

#endif
