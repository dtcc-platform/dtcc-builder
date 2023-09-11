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
  BoundingBox2D bounding_box{};

  /// Number of vertices (grid points) along X-axis
  size_t xsize{};

  /// Number of vertices (grid points) along Y-axis
  size_t ysize{};

  /// Resolution (grid size) along X-axis
  double xstep{};

  /// Resolution (grid size) along Y-axis
  double ystep{};

  /// Create empty grid
  Grid() = default;
  virtual ~Grid() {} // make the destructor virtual

  /// Create grid for given bounding box and size.
  ///
  /// @param bounding_box Bounding box for grid
  /// @param xsize Number of vertices (grid points) along X-axis
  /// @param y_size Number of vertices (grid points) along Y-axis
  Grid(const BoundingBox2D &bounding_box, size_t x_size, size_t y_size)
      : bounding_box(bounding_box), xsize(x_size), ysize(y_size)
  {
    assert(x_size > 1);
    assert(y_size > 1);
    xstep =
        (bounding_box.Q.x - bounding_box.P.x) / static_cast<double>(x_size - 1);
    ystep =
        (bounding_box.Q.y - bounding_box.P.y) / static_cast<double>(y_size - 1);
  }

  /// Compute total number of vertices
  ///
  /// @return Total number of vertices
  size_t num_vertices() const { return xsize * ysize; }

  /// Compute total number of cells
  ///
  /// @return Total number of cells
  size_t num_cells() const
  {
    if (xsize == 0 || ysize == 0)
      return 0;
    return (xsize - 1) * (ysize - 1);
  }

  /// Map vertex index to point.
  ///
  /// @param i Vertex index
  /// @return Vertex coordinates as a point
  Vector2D index_to_point(size_t i) const
  {
    const size_t ix = i % xsize;
    const size_t iy = i / xsize;
    return {bounding_box.P.x + ix * xstep, bounding_box.P.y + iy * ystep};
  }

  /// Map vertex index to (at most) 4 neighoring vertex indices.
  /// For efficiency, reserve return vector to size 4.
  ///
  /// @param i Vertex index
  /// @param indices Neighboring vertex indices
  void index_to_boundary(size_t i, std::vector<size_t> &indices) const
  {
    const size_t ix = i % xsize;
    const size_t iy = i / xsize;
    if (ix > 0)
      indices.push_back(i - 1);
    if (ix < xsize - 1)
      indices.push_back(i + 1);
    if (iy > 0)
      indices.push_back(i - xsize);
    if (iy < ysize - 1)
      indices.push_back(i + xsize);
  }

  /// Map vertex index to (at most) 4 neighoring vertex indices.
  ///
  /// @param i Vertex index
  /// @return Neighboring vertex indices
  std::vector<size_t> index_to_boundary(size_t i) const
  {
    std::vector<size_t> indices;
    indices.reserve(4);
    index_to_boundary(i, indices);
    return indices;
  }
  /// Map vertex index to (at most) 8 neighoring vertex indices.
  /// For efficiency, reserve return vector to size 4.
  ///
  /// @param i Vertex index
  /// @param indices Neighboring vertex indices
  void index_to_boundary8(size_t i, std::vector<size_t> &indices) const
  {
    const size_t ix = i % xsize;
    const size_t iy = i / xsize;
    if (ix > 0)
      indices.push_back(i - 1);
    if (ix < xsize - 1)
      indices.push_back(i + 1);
    if (iy > 0)
    {
      indices.push_back(i - xsize);
      const size_t ix2 = (i - xsize) % xsize;
      if (ix2 > 0)
        indices.push_back(i - xsize - 1);
      if (ix2 < xsize - 1)
        indices.push_back(i - xsize + 1);
    }
    if (iy < ysize - 1)
    {
      indices.push_back(i + xsize);
      const size_t ix2 = (i + xsize) % xsize;
      if (ix2 > 0)
        indices.push_back(i + xsize - 1);
      if (ix2 < xsize - 1)
        indices.push_back(i + xsize + 1);
    }
  }

  /// Map vertex index to (at most) 8 neighoring vertex indices.
  ///
  /// @param i Vertex index
  /// @return Neighboring vertex indices
  std::vector<size_t> index_to_boundary8(size_t i) const
  {
    std::vector<size_t> indices;
    indices.reserve(8);
    index_to_boundary(i, indices);
    return indices;
  }

  /// Map x and y indices to global index
  long int index_to_index(long int ix, long int iy) const
  {
    return ix + iy * xsize;
  }

  /// Map point to index of closest vertex.
  ///
  /// @param ix Grid index for x-direction (output)
  /// @param iy Grid index for y-direction (output)
  /// @param p Point
  void point_to_index(long int &ix, long int &iy, const Vector2D &p) const
  {
    const double _x = p.x - bounding_box.P.x;
    const double _y = p.y - bounding_box.P.y;
    ix = Utils::crop(std::lround(_x / xstep), xsize);
    iy = Utils::crop(std::lround(_y / ystep), ysize);
  }

  /// Map point to index of closest vertex.
  ///
  /// @param p Point
  /// @return Vertex index
  size_t point_to_index(const Vector2D &p) const
  {
    long int ix{}, iy{};
    point_to_index(ix, iy, p);
    return index_to_index(ix, iy);
  }

  /// Map point to cell and local coordinates.
  ///
  /// @param p Point
  /// @param i Index of cell containing point
  /// @param x Local X-coordinate on cell (scaled)
  /// @param y Local Y-coordinate on cell (scaled)
  void point_to_cell(const Vector2D &p, size_t &i, double &x, double &y) const
  {
    // Check that point is inside domain
    if (!Geometry::bounding_box_contains_2d(bounding_box, p))
      error("Point p = " + str(p) +
            " is outside of domain = " + str(bounding_box));

    // Compute grid cell containing point (lower left corner)
    const double _x = p.x - bounding_box.P.x;
    const double _y = p.y - bounding_box.P.y;
    const long int ix = Utils::crop(std::lround(_x / xstep - 0.5), xsize, 1);
    const long int iy = Utils::crop(std::lround(_y / ystep - 0.5), ysize, 1);
    i = iy * xsize + ix;

    // Map coordinates to [0, 1] x [0, 1] within grid cell
    x = (_x - ix * xstep) / xstep;
    y = (_y - iy * ystep) / ystep;
    assert(x >= 0.0 - Constants::epsilon);
    assert(x <= 1.0 + Constants::epsilon);
    assert(y >= 0.0 - Constants::epsilon);
    assert(y <= 1.0 + Constants::epsilon);
    if (x <= 0.0)
      x = 0.0;
    if (x >= 1.0)
      x = 1.0;
    if (y <= 0.0)
      y = 0.0;
    if (y >= 1.0)
      y = 1.0;
  }

  /// interpolate value on [0, 1] x [0, 1] using bilinear interpolation.
  ///
  /// @param x Local X-coordinate on [0, 1]
  /// @param y Local Y-coordinate on [0, 1]
  /// @param v00 Value at (0, 0)
  /// @param v10 Value at (1, 0)
  /// @param v01 Value at (0, 1)
  /// @param v11 Value at (0, 1)
  /// @return Value at (x, y)
  static double interpolate(
      double x, double y, double v00, double v10, double v01, double v11)
  {
    return (1.0 - x) * (1.0 - y) * v00 + x * (1.0 - y) * v10 +
           (1.0 - x) * y * v01 + x * y * v11;
  }

  /// Pretty-print
  std::string __str__() const override
  {
    return "Grid grid on " + str(bounding_box) + " of dimension " + str(xsize) +
           " x " + str(ysize);
  }
};
} // namespace DTCC_BUILDER

#endif
