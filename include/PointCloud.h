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

    /// Array of colors (one per point)
    std::vector<Color> Colors{};

    /// Array of classifications (one per point)
    std::vector<uint8_t> Classification{};

    /// Bounding box
    BoundingBox2D BoundingBox{};

    /// Create empty point cloud
    PointCloud() = default;

    /// Return density of point cloud (points per square meter)
    double Density() const
    {
      return static_cast<double>(Points.size()) / BoundingBox.Area();
    }

    /// set a default color to all points
    void InitColors(const Color &c)
    {
      Colors.clear();
      for (size_t i = 0; i < Points.size(); i++)
      {
        Colors.push_back(c);
      }
    }

    /// Set a default color to all points
    void InitClassification(uint8_t c)
    {
      Classification.clear();
      for (size_t i = 0; i < Points.size(); i++)
      {
        Classification.push_back(c);
      }

    }

    /// Clear all data
    void Clear()
    {
      Points.clear();
      Colors.clear();
      Classification.clear();
    }

    /// Pretty-print
    std::string __str__() const override
    {
      return "Point cloud on " + str(BoundingBox) + " with " +
             str(Points.size()) + " points and density " + str(Density()) +
             " m^-2";
    }
  };

} // namespace DTCC

#endif
