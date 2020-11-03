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
    PointCloud() = default;

    // set a default color to all points
    void InitColors(const Color &c)
    {
      for (size_t i = 0; i < Points.size(); i++)
      {
        Colors.push_back(c);
      }
    }
  
  void clear() 
    {
      Points.clear();
      Colors.clear();
    }

    /// Pretty-print
    std::string __str__() const override
    {
      return "Point cloud on " + str(BoundingBox) + " with " +
             str(Points.size()) + " points ";
    }
  };

} // namespace DTCC

#endif
