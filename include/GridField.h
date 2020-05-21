// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_GRID_FIELD_H
#define DTCC_GRID_FIELD_H

#include <assert.h>

#include "Grid.h"
#include "Field.h"
#include "Geometry.h"

namespace DTCC
{

  /// GridField2D represents a scalar field on a uniform 2D grid.
  /// The field can be efficiently evaluated at arbitrary points inside
  /// the grid domain. The value is computed by bilinear interpolation.
  class GridField2D : public Field2D
  {
  public:

    /// The grid
    Grid2D Grid{};

    /// Array of values (vertex values)
    std::vector<double> Values{};

    /// Create empty field
    GridField2D() {}

    /// Create zero field on given grid.
    ///
    /// @param grid The grid
    GridField2D(const Grid2D& grid)
      : Grid(grid)
    {
      // Initialize values to zero
      Values.resize(grid.NumVertices());
      std::fill(Values.begin(), Values.end(), 0.0);
    }

    /// Evaluate field at given point.
    ///
    /// @param p The point
    /// @return Value at point
    double operator()(const Point2D& p) const
    {
      // Check that point is inside domain
      if (!Geometry::BoundingBoxContains2D(Grid.BoundingBox, p))
      {
        std::cout << "p = " << p << std::endl;
        throw std::runtime_error("Unable to evaluate field; point outside of domain.");
      }

      // Compute grid cell containing point (lower left corner)
      const double _x = p.x - Grid.BoundingBox.P.x;
      const double _y = p.y - Grid.BoundingBox.P.y;
      const size_t ix = std::floor(_x / Grid.XStep);
      const size_t iy = std::floor(_y / Grid.YStep);
      const size_t i = iy*Grid.XSize + ix;
      assert(ix < Grid.XSize);
      assert(iy < Grid.YSize);
      assert(i < Values.size());

      // Map coordinates to [0, 1] x [0, 1] within grid cell
      const double X = (_x - ix*Grid.XStep) / Grid.XStep;
      const double Y = (_y - iy*Grid.YStep) / Grid.YStep;
      assert(X >= 0.0);
      assert(Y >= 0.0);
      assert(X <= 1.0);
      assert(Y <= 1.0);

      // Extract grid data
      const double z00 = Values[i];
      const double z10 = Values[i + 1];
      const double z01 = Values[i + Grid.XSize];
      const double z11 = Values[i + Grid.XSize + 1];

      // Compute value by bilinear interpolation
      return
        (1.0 - X) * (1.0 - Y) * z00 +
        (1.0 - X) * Y * z01 +
        X * (1.0 - Y) * z10 +
        X * Y * z11;
    }

  };

  /// GridField3D represents a scalar field on a uniform 3D grid.
  /// The field can be efficiently evaluated at arbitrary points inside
  /// the grid domain. The value is computed by trilinear interpolation.
  class GridField3D : public Field3D
  {
  public:

    /// The grid
    Grid3D Grid{};

    /// Array of values (vertex values)
    std::vector<double> Values{};

    /// Evaluate field at given point.
    ///
    /// @param p The point
    /// @return Value at point
    double operator()(const Point3D& p) const
    {
      // FIXME: In progress
      return 0.0;
    }

  };
}

#endif
