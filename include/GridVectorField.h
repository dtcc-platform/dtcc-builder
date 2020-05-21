// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_GRID_VECTOR_FIELD_H
#define DTCC_GRID_VECTOR_FIELD_H

#include "Grid.h"
#include "VectorField.h"
#include "Geometry.h"

namespace DTCC
{

  /// GridVectorField2D represents a vector field on a uniform 2D grid.
  /// The field can be efficiently evaluated at arbitrary points inside
  /// the grid domain. The value is computed by bilinear interpolation.
  class GridVectorField2D : public VectorField2D
  {
  public:

    /// The grid
    Grid2D Grid{};

    /// Array of values (flattened vertex values)
    std::vector<double> Values{};

    /// Evaluate field at given point.
    ///
    /// @param p The point
    /// @return Value at point
    Point2D operator()(const Point2D& p) const
    {
      // FIXME: In progress
      return Point2D();
    }

  };

  /// GridVectorField3D represents a vector field on a uniform 3D grid.
  /// The field can be efficiently evaluated at arbitrary points inside
  /// the grid domain. The value is computed by trilinear interpolation.
  class GridVectorField3D : public VectorField3D
  {
  public:

    /// The grid
    Grid3D Grid{};

    /// Array of values (flattened vertex values)
    std::vector<double> Values{};

    /// Evaluate field at given point.
    ///
    /// @param p The point
    /// @return Value at point
    Point3D operator()(const Point3D& p) const
    {
      // FIXME: In progress
      return Point3D();
    }

  };
}

#endif
