// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_POINT_CLOUD_H
#define DTCC_POINT_CLOUD_H

#include <vector>

#include "Color.h"
#include "Vector.h"
#include "BoundingBox.h"
#include "Logging.h"

namespace DTCC
{

  class PointCloud : public Printable
  {
  public:

    /// Array of points
    std::vector<Vector3D> Points{};

    std::vector<Color> Colors{};

    /// Bounding box
    BoundingBox2D BoundingBox{};

    // Create empty point cloud
  PointCloud() {}

    /// Pretty-print
    std::string __str__() const
    {
      return "Point cloud on " + str(BoundingBox) + " with " + str(Points.size()) + " points ";
    }

  };

} // namespace DTCC

#endif
