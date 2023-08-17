// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_GRID_H
#define DTCC_GRID_H

#include <cassert>

#include "BoundingBox.h"
#include "Geometry.h"
#include "Logging.h"
#include "Utils.h"

namespace DTCC_BUILDER
{

/// Grid represents a uniform 2D grid defined by a bounding box partitioned
/// into cells of equal size.
class Grid : public Printable
{
public:
  /// Bounding box for grid
  BoundingBox2D BoundingBox{};

  /// Number of vertices (grid points) along X-axis
  size_t XSize{};

  /// Number of vertices (grid points) along Y-axis
  size_t YSize{};

  /// Resolution (grid size) along X-axis
  double XStep{};

  /// Resolution (grid size) along Y-axis
  double YStep{};

  /// Create empty grid
  Grid() = default;

  /// Create grid for given bounding box and size.
  ///
  /// @param boundingBox Bounding box for grid
  /// @param xSize Number of vertices (grid points) along X-axis
  /// @param ySize Number of vertices (grid points) along Y-axis
  Grid(const BoundingBox2D &boundingBox, size_t xSize, size_t ySize)
      : BoundingBox(boundingBox), XSize(xSize), YSize(ySize)
  {
    assert(xSize > 1);
    assert(ySize > 1);
    XStep =
        (boundingBox.Q.x - boundingBox.P.x) / static_cast<double>(xSize - 1);
    YStep =
        (boundingBox.Q.y - boundingBox.P.y) / static_cast<double>(ySize - 1);
  }

  /// Compute total number of vertices
  ///
  /// @return Total number of vertices
  size_t NumVertices() const { return XSize * YSize; }

  /// Compute total number of cells
  ///
  /// @return Total number of cells
  size_t NumCells() const
  {
    if (XSize == 0 || YSize == 0)
      return 0;
    return (XSize - 1) * (YSize - 1);
  }

  /// Map vertex index to point.
  ///
  /// @param i Vertex index
  /// @return Vertex coordinates as a point
  Point2D Index2Point(size_t i) const
  {
    const size_t ix = i % XSize;
    const size_t iy = i / XSize;
    return {BoundingBox.P.x + ix * XStep, BoundingBox.P.y + iy * YStep};
  }

  /// Map vertex index to (at most) 4 neighoring vertex indices.
  /// For efficiency, reserve return vector to size 4.
  ///
  /// @param i Vertex index
  /// @param indices Neighboring vertex indices
  void Index2Boundary(size_t i, std::vector<size_t> &indices) const
  {
    const size_t ix = i % XSize;
    const size_t iy = i / XSize;
    if (ix > 0)
      indices.push_back(i - 1);
    if (ix < XSize - 1)
      indices.push_back(i + 1);
    if (iy > 0)
      indices.push_back(i - XSize);
    if (iy < YSize - 1)
      indices.push_back(i + XSize);
  }

  /// Map vertex index to (at most) 4 neighoring vertex indices.
  ///
  /// @param i Vertex index
  /// @return Neighboring vertex indices
  std::vector<size_t> Index2Boundary(size_t i) const
  {
    std::vector<size_t> indices;
    indices.reserve(4);
    Index2Boundary(i, indices);
    return indices;
  }
  /// Map vertex index to (at most) 8 neighoring vertex indices.
  /// For efficiency, reserve return vector to size 4.
  ///
  /// @param i Vertex index
  /// @param indices Neighboring vertex indices
  void Index2Boundary8(size_t i, std::vector<size_t> &indices) const
  {
    const size_t ix = i % XSize;
    const size_t iy = i / XSize;
    if (ix > 0)
      indices.push_back(i - 1);
    if (ix < XSize - 1)
      indices.push_back(i + 1);
    if (iy > 0)
    {
      indices.push_back(i - XSize);
      const size_t ix2 = (i - XSize) % XSize;
      if (ix2 > 0)
        indices.push_back(i - XSize - 1);
      if (ix2 < XSize - 1)
        indices.push_back(i - XSize + 1);
    }
    if (iy < YSize - 1)
    {
      indices.push_back(i + XSize);
      const size_t ix2 = (i + XSize) % XSize;
      if (ix2 > 0)
        indices.push_back(i + XSize - 1);
      if (ix2 < XSize - 1)
        indices.push_back(i + XSize + 1);
    }
  }

  /// Map vertex index to (at most) 8 neighoring vertex indices.
  ///
  /// @param i Vertex index
  /// @return Neighboring vertex indices
  std::vector<size_t> Index2Boundary8(size_t i) const
  {
    std::vector<size_t> indices;
    indices.reserve(8);
    Index2Boundary(i, indices);
    return indices;
  }

  /// Map x and y indices to global index
  long int Index2Index(long int ix, long int iy) const
  {
    return ix + iy * XSize;
  }

  /// Map point to index of closest vertex.
  ///
  /// @param ix Grid index for x-direction (output)
  /// @param iy Grid index for y-direction (output)
  /// @param p Point
  void Point2Index(long int &ix, long int &iy, const Point2D &p) const
  {
    const double _x = p.x - BoundingBox.P.x;
    const double _y = p.y - BoundingBox.P.y;
    ix = Utils::crop(std::lround(_x / XStep), XSize);
    iy = Utils::crop(std::lround(_y / YStep), YSize);
  }

  /// Map point to index of closest vertex.
  ///
  /// @param p Point
  /// @return Vertex index
  size_t Point2Index(const Point2D &p) const
  {
    long int ix{}, iy{};
    Point2Index(ix, iy, p);
    return Index2Index(ix, iy);
  }

  /// Map point to cell and local coordinates.
  ///
  /// @param p Point
  /// @param i Index of cell containing point
  /// @param x Local X-coordinate on cell (scaled)
  /// @param y Local Y-coordinate on cell (scaled)
  void Point2Cell(const Point2D &p, size_t &i, double &x, double &y) const
  {
    // Check that point is inside domain
    if (!Geometry::BoundingBoxContains2D(BoundingBox, p))
      error("Point p = " + str(p) +
            " is outside of domain = " + str(BoundingBox));

    // Compute grid cell containing point (lower left corner)
    const double _x = p.x - BoundingBox.P.x;
    const double _y = p.y - BoundingBox.P.y;
    const long int ix = Utils::crop(std::lround(_x / XStep - 0.5), XSize, 1);
    const long int iy = Utils::crop(std::lround(_y / YStep - 0.5), YSize, 1);
    i = iy * XSize + ix;

    // Map coordinates to [0, 1] x [0, 1] within grid cell
    x = (_x - ix * XStep) / XStep;
    y = (_y - iy * YStep) / YStep;
    assert(x >= 0.0 - Constants::Epsilon);
    assert(x <= 1.0 + Constants::Epsilon);
    assert(y >= 0.0 - Constants::Epsilon);
    assert(y <= 1.0 + Constants::Epsilon);
    if (x <= 0.0)
      x = 0.0;
    if (x >= 1.0)
      x = 1.0;
    if (y <= 0.0)
      y = 0.0;
    if (y >= 1.0)
      y = 1.0;
  }

  /// Interpolate value on [0, 1] x [0, 1] using bilinear interpolation.
  ///
  /// @param x Local X-coordinate on [0, 1]
  /// @param y Local Y-coordinate on [0, 1]
  /// @param v00 Value at (0, 0)
  /// @param v10 Value at (1, 0)
  /// @param v01 Value at (0, 1)
  /// @param v11 Value at (0, 1)
  /// @return Value at (x, y)
  static double Interpolate(
      double x, double y, double v00, double v10, double v01, double v11)
  {
    return (1.0 - x) * (1.0 - y) * v00 + x * (1.0 - y) * v10 +
           (1.0 - x) * y * v01 + x * y * v11;
  }

  /// Pretty-print
  std::string __str__() const override
  {
    return "Grid grid on " + str(BoundingBox) + " of dimension " + str(XSize) +
           " x " + str(YSize);
  }
};
} // namespace DTCC_BUILDER

#endif
