// Copyright (C) 2021 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_PLOTTING_H
#define DTCC_PLOTTING_H

#include "Logging.h"
#include "datamodel/Building.h"

namespace DTCCBUILDER
{
/// Simple utilities for printing Python plotting code for debugging
class Plotting
{
public:
  /// Print initialization code for plotting
  static void Init()
  {
    info("from pylab import *");
    info("from mpl_toolkits.mplot3d import Axes3D");
    info("ax = Axes3D(figure())");
  }

  /// Plot points (2D)
  static void Plot(const std::vector<Point2D> &points)
  {
    // Get coordinates
    const std::string x = GetX(points);
    const std::string y = GetY(points);

    // Print plot commands
    info("ax.scatter(" + x + ", " + y + ")");
  }

  /// Plot points (3D)
  static void Plot(const std::vector<Point3D> &points)
  {
    // Get coordinates
    const std::string x = GetX(points);
    const std::string y = GetY(points);
    const std::string z = GetZ(points);

    // Print plot commands
    info("ax.scatter(" + x + ", " + y + ", " + z + ")");
  }

  /// Plot building
  static void Plot(const Building &building, bool plotUUID = false)
  {
    // Build 3D points from footprint and ground height
    std::vector<Point3D> points;
    for (const auto &p2D : building.Footprint.Vertices)
    {
      Point3D p3D(p2D.x, p2D.y, building.GroundHeight);
      points.push_back(p3D);
    }

    // Get coordinates of footprint
    const std::string x = GetX(points, true);
    const std::string y = GetY(points, true);
    const std::string z = GetZ(points, true);

    // Get coordinates of ground points
    const std::string gx = GetX(building.GroundPoints);
    const std::string gy = GetY(building.GroundPoints);
    const std::string gz = GetZ(building.GroundPoints);

    // Get coordinates of roof points
    const std::string rx = GetX(building.RoofPoints);
    const std::string ry = GetY(building.RoofPoints);
    const std::string rz = GetZ(building.RoofPoints);

    // Get center
    Point2D c2D = Geometry::PolygonCenter2D(building.Footprint);
    Point3D c3D(c2D.x, c2D.y, building.GroundHeight);

    // Print plot commands
    info("ax.plot(" + x + ", " + y + ", " + z + ")");
    info("ax.scatter(" + gx + ", " + gy + ", " + gz + ", marker='o')");
    info("ax.scatter(" + rx + ", " + ry + ", " + rz + ", marker='.')");
    if (plotUUID)
      info("ax.text(" + str(c3D.x) + "," + str(c3D.y) + "," + str(c3D.z) +
           ",'" + building.UUID + "')");
  }

private:
  // FIXME: Would be nice to have one common function instead of three.

  // Get x-coordinates
  template <class T> static std::string GetX(T points, bool polygon = false)
  {
    size_t n = points.size();
    if (polygon)
      n += 1;
    std::string c = "[";
    for (size_t i = 0; i < n; i++)
    {
      c += str(points[i % points.size()].x);
      if (i < n - 1)
        c += ",";
    }
    c += "]";
    return c;
  }

  // Get y-coordinates
  template <class T> static std::string GetY(T points, bool polygon = false)
  {
    size_t n = points.size();
    if (polygon)
      n += 1;
    std::string c = "[";
    for (size_t i = 0; i < n; i++)
    {
      c += str(points[i % points.size()].y);
      if (i < n - 1)
        c += ",";
    }
    c += "]";
    return c;
  }

  // Get x-coordinates
  template <class T> static std::string GetZ(T points, bool polygon = false)
  {
    size_t n = points.size();
    if (polygon)
      n += 1;
    std::string c = "[";
    for (size_t i = 0; i < n; i++)
    {
      c += str(points[i % points.size()].z);
      if (i < n - 1)
        c += ",";
    }
    c += "]";
    return c;
  }
};
} // namespace DTCCBUILDER

#endif
