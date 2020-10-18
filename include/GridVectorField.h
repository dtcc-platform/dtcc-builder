// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_GRID_VECTOR_FIELD_H
#define DTCC_GRID_VECTOR_FIELD_H

#include "Grid.h"
#include "VectorField.h"
#include "Geometry.h"
#include "Logging.h"

namespace DTCC
{

  /// GridVectorField2D represents a vector field on a uniform 2D grid.
  /// The field can be efficiently evaluated at arbitrary points inside
  /// the grid domain. The value is computed by bilinear interpolation.
  class GridVectorField2D : public VectorField2D, public Printable
  {
  public:

    /// The grid
    Grid2D Grid{};

    /// Array of values (flattened vertex values: [x y x y ...])
    std::vector<double> Values{};

    /// Create empty field
    GridVectorField2D() = default;

    /// Create zero field on given grid.
    ///
    /// @param grid The grid
    explicit GridVectorField2D(const Grid2D& grid) : Grid(grid)
    {
      // Initialize values to zero
      Values.resize(2*grid.NumVertices());
      std::fill(Values.begin(), Values.end(), 0.0);
    }

    /// Evaluate field at given point.
    ///
    /// @param p The point
    /// @return Value at point
    Vector2D operator()(const Point2D& p) const override
    {
      // Map point to cell
      size_t i{};
      double x{}, y{};
      Grid.Point2Cell(p, i, x, y);

      // Extract grid data
      const double v00x = Values[2*i];
      const double v10x = Values[2*(i + 1)];
      const double v01x = Values[2*(i + Grid.XSize)];
      const double v11x = Values[2*(i + Grid.XSize + 1)];
      const double v00y = Values[2*i + 1];
      const double v10y = Values[2*(i + 1) + 1];
      const double v01y = Values[2*(i + Grid.XSize) + 1];
      const double v11y = Values[2*(i + Grid.XSize + 1) + 1];

      // Compute value by bilinear interpolation
      Vector2D v{};
      v.x = Grid.Interpolate(x, y, v00x, v10x, v01x, v11x);
      v.y = Grid.Interpolate(x, y, v00y, v10y, v01y, v11y);

      return v;
    }

    /// Interpolate given field at vertices.
    ///
    /// @param field The field to be interpolated
    void Interpolate(const VectorField2D& field)
    {
      // Iterate over vertices and evaluate field
      for (size_t i = 0; i < Grid.NumVertices(); i++)
      {
        const Vector2D v = field(Grid.Index2Point(i));
        Values[2*i] = v.x;
        Values[2*i + 1] = v.y;
      }
    }

    /// Pretty-print
    std::string __str__() const override
    {
      return "2D vector field on " + str(Grid);
    }

  };

  /// GridVectorField3D represents a vector field on a uniform 3D grid.
  /// The field can be efficiently evaluated at arbitrary points inside
  /// the grid domain. The value is computed by trilinear interpolation.
  class GridVectorField3D : public VectorField3D, public Printable
  {
  public:

    /// The grid
    Grid3D Grid{};

    /// Array of values (flattened vertex values: [x y z x y z ...])
    std::vector<double> Values{};

    /// Create empty field
    GridVectorField3D() = default;

    /// Create zero field on given grid.
    ///
    /// @param grid The grid
    explicit GridVectorField3D(const Grid3D& grid) : Grid(grid)
    {
      // Initialize values to zero
      Values.resize(3*grid.NumVertices());
      std::fill(Values.begin(), Values.end(), 0.0);
    }

    /// Evaluate field at given point.
    ///
    /// @param p The point
    /// @return Value at point
    Vector3D operator()(const Point3D& p) const override
    {
      // Map point to cell
      size_t i{};
      double x{}, y{}, z{};
      Grid.Point2Cell(p, i, x, y, z);

      // Extract grid data
      const double v000x = Values[3*i];
      const double v100x = Values[3*(i + 1)];
      const double v010x = Values[3*(i + Grid.XSize)];
      const double v110x = Values[3*(i + Grid.XSize + 1)];
      const double v001x = Values[3*(i + Grid.XSize * Grid.YSize)];
      const double v101x = Values[3*(i + 1 + Grid.XSize * Grid.YSize)];
      const double v011x = Values[3*(i + Grid.XSize + Grid.XSize * Grid.YSize)];
      const double v111x = Values[3*(i + Grid.XSize + 1 + Grid.XSize * Grid.YSize)];
      const double v000y = Values[3*i + 1];
      const double v100y = Values[3*(i + 1) + 1];
      const double v010y = Values[3*(i + Grid.XSize) + 1];
      const double v110y = Values[3*(i + Grid.XSize + 1) + 1];
      const double v001y = Values[3*(i + Grid.XSize * Grid.YSize) + 1];
      const double v101y = Values[3*(i + 1 + Grid.XSize * Grid.YSize) + 1];
      const double v011y = Values[3*(i + Grid.XSize + Grid.XSize * Grid.YSize) + 1];
      const double v111y = Values[3*(i + Grid.XSize + 1 + Grid.XSize * Grid.YSize) + 1];
      const double v000z = Values[3*i + 2];
      const double v100z = Values[3*(i + 1) + 2];
      const double v010z = Values[3*(i + Grid.XSize) + 2];
      const double v110z = Values[3*(i + Grid.XSize + 1) + 2];
      const double v001z = Values[3*(i + Grid.XSize * Grid.YSize) + 2];
      const double v101z = Values[3*(i + 1 + Grid.XSize * Grid.YSize) + 2];
      const double v011z = Values[3*(i + Grid.XSize + Grid.XSize * Grid.YSize) + 2];
      const double v111z = Values[3*(i + Grid.XSize + 1 + Grid.XSize * Grid.YSize) + 2];

      // Compute value by trilinear interpolation
      Vector3D v;
      v.x = Grid.Interpolate(x, y, z, v000x, v100x, v010x, v110x, v001x, v101x, v011x, v111x);
      v.y = Grid.Interpolate(x, y, z, v000y, v100y, v010y, v110y, v001y, v101y, v011y, v111y);
      v.z = Grid.Interpolate(x, y, z, v000z, v100z, v010z, v110z, v001z, v101z, v011z, v111z);

      return v;
    }

    /// Interpolate given field at vertices.
    ///
    /// @param field The field to be interpolated
    void Interpolate(const VectorField3D& field)
    {
      // Iterate over vertices and evaluate field
      for (size_t i = 0; i < Grid.NumVertices(); i++)
      {
        const Vector3D v = field(Grid.Index2Point(i));
        Values[3*i] = v.x;
        Values[3*i + 1] = v.y;
        Values[3*i + 2] = v.z;
      }
    }

    /// Pretty-print
    std::string __str__() const override
    {
      return "3D vector field on " + str(Grid);
    }

  };
}

#endif
