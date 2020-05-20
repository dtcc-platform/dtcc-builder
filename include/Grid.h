// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_GRID_H
#define DTCC_GRID_H

#include <assert.h>

#include "BoundingBox.h"

namespace DTCC
{

  /// Grid2D represents a uniform 2D grid defined by a bounding box partitioned
  /// into cells of equal size.
  class Grid2D
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
    Grid2D() {}

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

  };

  /// Grid3D represents a uniform 3D grid defined by a bounding box partitioned
  /// into cells of equal size.
  class Grid3D
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
    Grid3D() {}

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

  };

}

#endif
