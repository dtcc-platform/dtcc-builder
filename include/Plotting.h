// Copyright (C) 2021 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_PLOTTING_H
#define DTCC_PLOTTING_H

#include "Building.h"
#include "Logging.h"

namespace DTCC
{
/// Simple utilities for printing Python plotting code for debugging
class Plotting
{
public:
  /// Print initialization code for plotting
  static void Init()
  {
    Info("from pylab import *");
    Info("from mpl_toolkits.mplot3d import Axes3D");
    Info("ax = Axes3D(figure())");
  }

  /// Print plotting code for points (2D)
  static void Plot(const std::vector<Point2D> points, std::string marker = ".r")
  {
    // Extract coordinates
    const size_t n = points.size();
    std::string x, y;
    for (size_t i = 0; i < n; i++)
    {
      if (i == 0)
      {
        x += "[";
        y += "[";
      }

      const Point2D &p = points[i];
      x += str(p.x);
      y += str(p.y);

      if (i == n - 1)
      {
        x += "]";
        y += "]";
      }
      else
      {
        x += ",";
        y += ",";
      }
    }

    // Print plot command
    Info("ax.plot(" + x + ", " + y + ", '" + marker + "')");
  }

  /// Print plotting code for points (3D)
  static void Plot(const std::vector<Point3D> points, std::string marker = ".r")
  {
    // Extract coordinates
    const size_t n = points.size();
    std::string x, y, z;
    for (size_t i = 0; i < n; i++)
    {
      if (i == 0)
      {
        x += "[";
        y += "[";
        z += "[";
      }

      const Point3D &p = points[i];
      x += str(p.x);
      y += str(p.y);
      z += str(p.z);

      if (i == n - 1)
      {
        x += "]";
        y += "]";
        z += "]";
      }
      else
      {
        x += ",";
        y += ",";
        z += ",";
      }
    }

    // Print plot command
    Info("ax.scatter(" + x + ", " + y + ", " + z + ", '" + marker + "')");
  }

  /// Print plotting code for building
  static void Plot(const Building &building)
  {
    Plot(building.Footprint.Vertices, "-bo");
  }
};
}

#endif
