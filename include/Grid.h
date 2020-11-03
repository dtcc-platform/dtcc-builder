// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_GRID_H
#define DTCC_GRID_H

#include <cassert>

#include "BoundingBox.h"
#include "Geometry.h"
#include "Utils.h"
#include "Logging.h"

namespace DTCC
{

  /// Grid2D represents a uniform 2D grid defined by a bounding box partitioned
  /// into cells of equal size.
  class Grid2D : public Printable
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
    Grid2D() = default;

    /// Create grid for given bounding box and size.
    ///
    /// @param boundingBox Bounding box for grid
    /// @param xSize Number of vertices (grid points) along X-axis
    /// @param ySize Number of vertices (grid points) along Y-axis
    Grid2D(const BoundingBox2D& boundingBox, size_t xSize, size_t ySize)
      : BoundingBox(boundingBox), XSize(xSize), YSize(ySize)
    {
      assert(xSize > 1);
      assert(ySize > 1);
      XStep = (boundingBox.Q.x - boundingBox.P.x) / static_cast<double>(xSize - 1);
      YStep = (boundingBox.Q.y - boundingBox.P.y) / static_cast<double>(ySize - 1);
    }

    /// Compute total number of vertices
    ///
    /// @return Total number of vertices
    size_t NumVertices() const
    {
      return XSize*YSize;
    }

    /// Compute total number of cells
    ///
    /// @return Total number of cells
    size_t NumCells() const
    {
      if (XSize == 0 || YSize == 0)
        return 0;
      return (XSize - 1)*(YSize - 1);
    }

    /// Map vertex index to point.
    ///
    /// @param i Vertex index
    /// @return Vertex coordinates as a point
    Point2D Index2Point(size_t i) const
    {
      const size_t ix = i % XSize;
      const size_t iy = i / XSize;
      return {BoundingBox.P.x + ix * XStep,
                     BoundingBox.P.y + iy * YStep};
    }

    /// Map vertex index to (at most) 4 neighoring vertex indices.
    /// For efficiency, reserve return vector to size 4.
    ///
    /// @param i Vertex index
    /// @param indices Neighboring vertex indices
    void Index2Boundary(size_t i, std::vector<size_t>& indices) const
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

    /// Map point to index of closest vertex.
    ///
    /// @param p Point
    /// @return Vertex index
    size_t Point2Index(const Point2D& p) const
    {
      const double _x = p.x - BoundingBox.P.x;
      const double _y = p.y - BoundingBox.P.y;
      const long int ix = Utils::crop(std::lround(_x / XStep), XSize);
      const long int iy = Utils::crop(std::lround(_y / YStep), YSize);
      return ix + iy * XSize;
    }

    /// Map point to cell and local coordinates.
    ///
    /// @param p Point
    /// @param i Index of cell containing point
    /// @param x Local X-coordinate on cell (scaled)
    /// @param y Local Y-coordinate on cell (scaled)
    void Point2Cell(const Point2D& p, size_t& i, double& x, double& y) const
    {
      // Check that point is inside domain
      if (!Geometry::BoundingBoxContains2D(BoundingBox, p))
        Error("Point p = " + str(p) + " is outside of domain = " + str(BoundingBox));

      // Compute grid cell containing point (lower left corner)
      const double _x = p.x - BoundingBox.P.x;
      const double _y = p.y - BoundingBox.P.y;
      const long int ix = Utils::crop(std::lround(_x / XStep - 0.5), XSize, 1);
      const long int iy = Utils::crop(std::lround(_y / YStep - 0.5), YSize, 1);
      i = iy * XSize + ix;

      // Map coordinates to [0, 1] x [0, 1] within grid cell
      x = (_x - ix*XStep) / XStep;
      y = (_y - iy*YStep) / YStep;
      assert(x >= 0.0);
      assert(y >= 0.0);
      assert(x <= 1.0);
      assert(y <= 1.0);
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
    static double Interpolate(double x, double y,
                              double v00, double v10, double v01, double v11)
    {
      return
        (1.0 - x) * (1.0 - y) * v00 +
        x * (1.0 - y) * v10 +
        (1.0 - x) * y * v01 +
        x * y * v11;
    }

    /// Pretty-print
    std::string __str__() const override
    {
      return "2D grid on " + str(BoundingBox) + " of dimension " +
        str(XSize) + " x " + str(YSize);
    }

  };

  /// Grid3D represents a uniform 3D grid defined by a bounding box partitioned
  /// into cells of equal size.
  class Grid3D : public Printable
  {
  public:

    /// Bounding box for grid
    BoundingBox3D BoundingBox{};

    /// Number of vertices (grid points) along X-axis
    size_t XSize{};

    /// Number of vertices (grid points) along Y-axis
    size_t YSize{};

    /// Number of vertices (grid points) along Z-axis
    size_t ZSize{};

    /// Resolution (grid size) along X-axis
    double XStep{};

    /// Resolution (grid size) along Y-axis
    double YStep{};

    /// Resolution (grid size) along Z-axis
    double ZStep{};

    /// Create empty grid
    Grid3D() = default;

    /// Create grid for given bounding box and size.
    ///
    /// @param boundingBox Bounding box for grid
    /// @param xSize Number of vertices (grid points) along X-axis
    /// @param ySize Number of vertices (grid points) along Y-axis
    /// @param zSize Number of vertices (grid points) along Z-axis
    Grid3D(const BoundingBox3D& boundingBox, size_t xSize, size_t ySize, size_t zSize)
      : BoundingBox(boundingBox), XSize(xSize), YSize(ySize), ZSize(zSize)
    {
      assert(xSize > 1);
      assert(ySize > 1);
      assert(zSize > 1);
      XStep = (boundingBox.Q.x - boundingBox.P.x) / static_cast<double>(xSize - 1);
      YStep = (boundingBox.Q.y - boundingBox.P.y) / static_cast<double>(ySize - 1);
      ZStep = (boundingBox.Q.z - boundingBox.P.z) / static_cast<double>(zSize - 1);
    }

    /// Compute total number of vertices
    ///
    /// @return Total number of vertices
    size_t NumVertices() const
    {
      return XSize*YSize*ZSize;
    }

    /// Compute total number of cells
    ///
    /// @return Total number of cells
    size_t NumCells() const
    {
      if (XSize == 0 || YSize == 0 || ZSize == 0)
        return 0;
      return (XSize - 1)*(YSize - 1)*(ZSize - 1);
    }

    /// Map vertex index to point.
    ///
    /// @param i Vertex index
    /// @return Vertex coordinates as a point
    Point3D Index2Point(size_t i) const
    {
      const size_t ix = i % XSize;
      const size_t iy = (i / XSize) % YSize;
      const size_t iz = i / (XSize * YSize);
      return {BoundingBox.P.x + ix * XStep,
                     BoundingBox.P.y + iy * YStep,
                     BoundingBox.P.z + iz * ZStep};
    }

    /// Map point to index of closest vertex.
    ///
    /// @param p Point
    /// @return Vertex index
    size_t Point2Index(const Point3D& p) const
    {
      const double _x = p.x - BoundingBox.P.x;
      const double _y = p.y - BoundingBox.P.y;
      const double _z = p.z - BoundingBox.P.z;
      const long int ix = Utils::crop(std::lround(_x / XStep), XSize);
      const long int iy = Utils::crop(std::lround(_y / YStep), YSize);
      const long int iz = Utils::crop(std::lround(_z / ZStep), ZSize);
      return ix + iy * XSize + iz * XSize * YSize;
    }

    /// Map point to cell and local coordinates.
    ///
    /// @param p Point
    /// @param i Index of cell containing point
    /// @param X Local X-coordinate on cell (scaled)
    /// @param Y Local Y-coordinate on cell (scaled)
    /// @param Z Local Z-coordinate on cell (scaled)
    void Point2Cell(const Point3D& p, size_t& i, double& x, double& y, double& z) const
    {
      // Check that point is inside domain
      if (!Geometry::BoundingBoxContains3D(BoundingBox, p))
        Error("Point p = " + str(p) + " is outside of domain = " + str(BoundingBox));

      // Compute grid cell containing point (lower left corner)
      const double _x = p.x - BoundingBox.P.x;
      const double _y = p.y - BoundingBox.P.y;
      const double _z = p.z - BoundingBox.P.z;
      const long int ix = Utils::crop(std::lround(_x / XStep - 0.5), XSize, 1);
      const long int iy = Utils::crop(std::lround(_y / YStep - 0.5), YSize, 1);
      const long int iz = Utils::crop(std::lround(_z / ZStep - 0.5), ZSize, 1);
      i = iz * XSize * YSize + iy * XSize + ix;

      // Map coordinates to [0, 1] x [0, 1] within grid cell
      x = (_x - ix * XStep) / XStep;
      y = (_y - iy * YStep) / YStep;
      z = (_z - iz * ZStep) / ZStep;
      assert(x >= 0.0);
      assert(y >= 0.0);
      assert(z >= 0.0);
      assert(x <= 1.0);
      assert(y <= 1.0);
      assert(z <= 1.0);
    }

    /// Interpolate value on [0, 1] x [0, 1] x [0, 1] using trilinear interpolation.
    ///
    /// @param x Local X-coordinate on [0, 1]
    /// @param y Local Y-coordinate on [0, 1]
    /// @param z Local Z-coordinate on [0, 1]
    /// @param v000 Value at (0, 0, 0)
    /// @param v100 Value at (1, 0, 0)
    /// @param v010 Value at (0, 1, 0)
    /// @param v110 Value at (0, 1, 0)
    /// @param v001 Value at (0, 0, 1)
    /// @param v101 Value at (1, 0, 1)
    /// @param v011 Value at (0, 1, 1)
    /// @param v111 Value at (0, 1, 1)
    /// @return Value at (x, y, z)
    static double Interpolate(double x, double y, double z,
                              double v000, double v100, double v010, double v110,
                              double v001, double v101, double v011, double v111)
    {
      return
        (1.0 - x) * (1.0 - y) * (1.0 - z) * v000 +
        x * (1.0 - y) * (1.0 - z) * v100 +
        (1.0 - x) * y * (1.0 - z) * v010 +
        x * y *  (1.0 - z) * v110 +
        (1.0 - x) * (1.0 - y) * z * v001 +
        x * (1.0 - y) * z * v101 +
        (1.0 - x) * y * z * v011 +
        x * y * z * v111;
    }

    /// Pretty-print
    std::string __str__() const override
    {
      return "3D grid on " + str(BoundingBox) + " of dimension " +
        str(XSize) + " x " + str(YSize) + " x " + str(YSize);
    }

  };

}

#endif
