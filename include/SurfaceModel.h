// Digital Surface Model (DSM) for elevation (including buildings).
// Copyright (C) 2019 Anders Logg.

#ifndef DTCC_SURFACE_MODEL_H
#define DTCC_SURFACE_MODEL_H

#include <vector>

#include "GeoReference.h"
#include "Vector.h"

namespace DTCC
{

class SurfaceModel
{
public:
  // Grid dimensions
  double XMin, XMax, YMin, YMax;

  // Resolution
  double Resolution;

  // Number of grid points
  size_t SizeX, SizeY;

  // Grid data (flattened array of (x, y) coordinates)
  std::vector<double> GridData;

  // Create empty height map
  SurfaceModel(
      double xMin, double xMax, double yMin, double yMax, double resolution)
      : XMin(xMin), XMax(xMax), YMin(yMin), YMax(yMax), Resolution(resolution)
  {
    // Initialize grid data
    SizeX = (XMax - XMin) / resolution + 1;
    SizeY = (YMax - YMin) / resolution + 1;
    GridData.resize(SizeX * SizeY);
  }

  // Return height (z) at 2D point p
  double operator()(const Vector2D &p) const { return (*this)(p.x, p.y); }

  // Return height (z) at 2D point (x, y)
  double operator()(double x, double y) const { return 0.0; }

  // Return coordinate number i (of flattened grid data)
  Vector2D Coordinate(size_t i) const
  {
    const size_t ix = i % SizeX;
    const size_t iy = i / SizeX;
    return Vector2D(ix * Resolution, iy * Resolution);
  }
};

} // namespace DTCC

#endif
