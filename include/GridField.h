// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_GRID_FIELD_H
#define DTCC_GRID_FIELD_H

#include <cassert>

#include "Grid.h"
#include "Field.h"
#include "Geometry.h"
#include "Logging.h"

namespace DTCC_BUILDER
{

  /// GridField2D represents a scalar field on a uniform 2D grid.
  /// The field can be efficiently evaluated at arbitrary points inside
  /// the grid domain. The value is computed by bilinear interpolation.
  class GridField2D : public Field2D, public Printable
  {
  public:

    /// The grid
    Grid2D Grid{};

    /// Array of values (vertex values)
    std::vector<double> Values{};

    /// Create empty field
    GridField2D() = default;

    /// Create zero field on given grid.
    ///
    /// @param grid The grid
    explicit GridField2D(const Grid2D& grid) : Grid(grid)
    {
      // Initialize values to zero
      Values.resize(grid.NumVertices());
      std::fill(Values.begin(), Values.end(), 0.0);
    }

    /// Evaluate field at given point.
    ///
    /// @param p The point
    /// @return Value at point
    double operator()(const Point2D& p) const override
    {
      // Map point to cell
      size_t i{};
      double x{}, y{};
      Grid.Point2Cell(p, i, x, y);

      // Extract grid data
      const double v00 = Values[i];
      const double v10 = Values[i + 1];
      const double v01 = Values[i + Grid.XSize];
      const double v11 = Values[i + Grid.XSize + 1];

      // Compute value by bilinear interpolation
      return Grid2D::Interpolate(x, y, v00, v10, v01, v11);
    }

    double Nearest(const Point2D& p) const
    {
      size_t i{};
      double x{}, y{};
      Grid.Point2Cell(p, i, x, y);
      return Values[i];
    }

    /// Interpolate given field at vertices.
    ///
    /// @param field The field to be interpolated
    void Interpolate(const Field2D& field)
    {
      // Iterate over vertices and evaluate field
      for (size_t i = 0; i < Grid.NumVertices(); i++)
        Values[i] = field(Grid.Index2Point(i));
    }

    /// Compute minimal vertex value.
    ///
    /// @return Minimal vertex value
    double Min() const
    {
      return *std::min_element(Values.begin(), Values.end());
    }

    /// Compute maximal of vertex value.
    ///
    /// @return Maximal vertex value
    double Max() const
    {
      return *std::max_element(Values.begin(), Values.end());
    }

    /// Compute mean vertex value
    ///
    /// @return Mean vertex value
    double Mean() const
    {
      double mean = 0.0;
      for (const auto &value : Values)
        mean += value;
      return mean / static_cast<double>(Values.size());
    }

    /// Pretty-print
    std::string __str__() const override
    {
      return "2D field on " + str(Grid);
    }

  };

  /// GridField3D represents a scalar field on a uniform 3D grid.
  /// The field can be efficiently evaluated at arbitrary points inside
  /// the grid domain. The value is computed by trilinear interpolation.
  class GridField3D : public Field3D, public Printable
  {
  public:

    /// The grid
    Grid3D Grid{};

    /// Array of values (vertex values)
    std::vector<double> Values{};

    /// Create empty field
    GridField3D() = default;

    /// Create zero field on given grid.
    ///
    /// @param grid The grid
    explicit GridField3D(const Grid3D& grid) : Grid(grid)
    {
      // Initialize values to zero
      Values.resize(grid.NumVertices());
      std::fill(Values.begin(), Values.end(), 0.0);
    }

    /// Evaluate field at given point.
    ///
    /// @param p The point
    /// @return Value at point
    double operator()(const Point3D& p) const override
    {
      // Map point to cell
      size_t i{};
      double x{}, y{}, z{};
      Grid.Point2Cell(p, i, x, y, z);

      // Extract grid data
      const double v000 = Values[i];
      const double v100 = Values[i + 1];
      const double v010 = Values[i + Grid.XSize];
      const double v110 = Values[i + Grid.XSize + 1];
      const double v001 = Values[i + Grid.XSize * Grid.YSize];
      const double v101 = Values[i + 1 + Grid.XSize * Grid.YSize];
      const double v011 = Values[i + Grid.XSize + Grid.XSize * Grid.YSize];
      const double v111 = Values[i + Grid.XSize + 1 + Grid.XSize * Grid.YSize];

      // Compute value by trilinear interpolation
      return Grid3D::Interpolate(x, y, z, v000, v100, v010, v110, v001, v101,
                                 v011, v111);
    }

    double Nearest(const Point3D& p) const 
    {
      size_t i{};
      double x{}, y{}, z{};
      Grid.Point2Cell(p, i, x, y, z);
      return  Values[i];
    }

    /// Interpolate given field at vertices.
    ///
    /// @param field The field to be interpolated
    void Interpolate(const Field3D& field)
    {
      // Iterate over vertices and evaluate field
      for (size_t i = 0; i < Grid.NumVertices(); i++)
        Values[i] = field(Grid.Index2Point(i));
    }

    /// Compute minimal vertex value.
    ///
    /// @return Minimal vertex value
    double Min() const
    {
      return *std::min_element(Values.begin(), Values.end());
    }

    /// Compute maximal of vertex value.
    ///
    /// @return Maximal vertex value
    double Max() const
    {
      return *std::max_element(Values.begin(), Values.end());
    }

    /// Compute mean vertex value
    ///
    /// @return Mean vertex value
    double Mean() const
    {
      double mean = 0.0;
      for (const auto &value : Values)
        mean += value;
      return mean / static_cast<double>(Values.size());
    }

    /// Pretty-print
    std::string __str__() const override
    {
      return "3D field on " + str(Grid);
    }

  };
}

#endif
