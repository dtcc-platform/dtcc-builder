// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_GRID_FIELD_H
#define DTCC_GRID_FIELD_H

#include <assert.h>

#include "Grid.h"
#include "Field.h"
#include "Geometry.h"
#include "Utils.h"
#include "Logging.h"

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
        Error("Point p = " + str(p) + " is outside of domain = " + str(Grid.BoundingBox));

      // Compute grid cell containing point (lower left corner)
      const double _x = p.x - Grid.BoundingBox.P.x;
      const double _y = p.y - Grid.BoundingBox.P.y;
      const long int ix = Utils::crop(std::lround(_x / Grid.XStep - 0.5), Grid.XSize, 1);
      const long int iy = Utils::crop(std::lround(_y / Grid.YStep - 0.5), Grid.YSize, 1);
      const size_t i = iy * Grid.XSize + ix;

      // Map coordinates to [0, 1] x [0, 1] within grid cell
      const double X = (_x - ix*Grid.XStep) / Grid.XStep;
      const double Y = (_y - iy*Grid.YStep) / Grid.YStep;
      assert(X >= 0.0);
      assert(Y >= 0.0);
      assert(X <= 1.0);
      assert(Y <= 1.0);

      // Extract grid data
      const double v00 = Values[i];
      const double v10 = Values[i + 1];
      const double v01 = Values[i + Grid.XSize];
      const double v11 = Values[i + Grid.XSize + 1];

      // Compute value by bilinear interpolation
      return
        (1.0 - X) * (1.0 - Y) * v00 +
        X * (1.0 - Y) * v10 +
        (1.0 - X) * Y * v01 +
        X * Y * v11;
    }

    /// Interpolate given field at vertices.
    ///
    /// @param field The field to be interpolated
    void Interpolate(const Field2D& field)
    {
      // Iterate over vertices and evaluate field
      for (size_t i = 0; i < Values.size(); i++)
        Values[i] = field(Index2Point(i));
    }

    /// Map vertex index to point.
    ///
    /// @param i Vertex index
    /// @return Vertex coordinates as a point
    Point2D Index2Point(size_t i) const
    {
      const size_t ix = i % Grid.XSize;
      const size_t iy = i / Grid.XSize;
      return Point2D(Grid.BoundingBox.P.x + ix * Grid.XStep,
                     Grid.BoundingBox.P.y + iy * Grid.YStep);
    }

    /// Map point to index of closest vertex.
    ///
    /// @param p Point
    /// @return Vertex index
    size_t Point2Index(const Point2D& p)
    {
      const double _x = p.x - Grid.BoundingBox.P.x;
      const double _y = p.y - Grid.BoundingBox.P.y;
      const long int ix = Utils::crop(std::lround(_x / Grid.XStep), Grid.XSize);
      const long int iy = Utils::crop(std::lround(_y / Grid.YStep), Grid.YSize);
      return ix + iy * Grid.XSize;
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

    /// Create zero field on given grid.
    ///
    /// @param grid The grid
    GridField3D(const Grid3D& grid)
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
    double operator()(const Point3D& p) const
    {
      // Check that point is inside domain
      if (!Geometry::BoundingBoxContains3D(Grid.BoundingBox, p))
        Error("Point p = " + str(p) + " is outside of domain = " + str(Grid.BoundingBox));

      // Compute grid cell containing point (lower left corner)
      const double _x = p.x - Grid.BoundingBox.P.x;
      const double _y = p.y - Grid.BoundingBox.P.y;
      const double _z = p.z - Grid.BoundingBox.P.z;
      const long int ix = Utils::crop(std::lround(_x / Grid.XStep - 0.5), Grid.XSize, 1);
      const long int iy = Utils::crop(std::lround(_y / Grid.YStep - 0.5), Grid.YSize, 1);
      const long int iz = Utils::crop(std::lround(_z / Grid.ZStep - 0.5), Grid.ZSize, 1);
      const size_t i = iz * Grid.XSize * Grid.YSize + iy * Grid.XSize + ix;

      // Map coordinates to [0, 1] x [0, 1] within grid cell
      const double X = (_x - ix * Grid.XStep) / Grid.XStep;
      const double Y = (_y - iy * Grid.YStep) / Grid.YStep;
      const double Z = (_z - iz * Grid.ZStep) / Grid.ZStep;
      assert(X >= 0.0);
      assert(Y >= 0.0);
      assert(Z >= 0.0);
      assert(X <= 1.0);
      assert(Y <= 1.0);
      assert(Z <= 1.0);

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
      return
        (1.0 - X) * (1.0 - Y) * (1.0 - Z) * v000 +
        X * (1.0 - Y) * (1.0 - Z) * v100 +
        (1.0 - X) * Y * (1.0 - Z) * v010 +
        X * Y *  (1.0 - Z) * v110 +
        (1.0 - X) * (1.0 - Y) * Z * v001 +
        X * (1.0 - Y) * Z * v101 +
        (1.0 - X) * Y * Z * v011 +
        X * Y * Z * v111;
    }

    /// Interpolate given field at vertices.
    ///
    /// @param field The field to be interpolated
    void Interpolate(const Field3D& field)
    {
      // Iterate over vertices and evaluate field
      for (size_t i = 0; i < Values.size(); i++)
        Values[i] = field(Index2Point(i));
    }

    /// Map vertex index to point.
    ///
    /// @param i Vertex index
    /// @return Vertex coordinates as a point
    Point3D Index2Point(size_t i) const
    {
      const size_t ix = i % Grid.XSize;
      const size_t iy = (i / Grid.XSize) % Grid.YSize;
      const size_t iz = i / (Grid.XSize * Grid.YSize);
      return Point3D(Grid.BoundingBox.P.x + ix * Grid.XStep,
                     Grid.BoundingBox.P.y + iy * Grid.YStep,
                     Grid.BoundingBox.P.z + iz * Grid.ZStep);
    }

    /// Map point to index of closest vertex.
    ///
    /// @param p Point
    /// @return Vertex index
    size_t Point2Index(const Point3D& p)
    {
      const double _x = p.x - Grid.BoundingBox.P.x;
      const double _y = p.y - Grid.BoundingBox.P.y;
      const double _z = p.z - Grid.BoundingBox.P.z;
      const long int ix = Utils::crop(std::lround(_x / Grid.XStep), Grid.XSize);
      const long int iy = Utils::crop(std::lround(_y / Grid.YStep), Grid.YSize);
      const long int iz = Utils::crop(std::lround(_z / Grid.ZStep), Grid.ZSize);
      return ix + iy * Grid.XSize + iz * Grid.XSize * Grid.YSize;
    }
  };
}

#endif
