// Point cloud (array of 3D points).
// Copyright (C) 2019 Anders Logg.

#ifndef VC_POINT_CLOUD_H
#define VC_POINT_CLOUD_H

#include <vector>

#include "Point.h"

namespace VirtualCity
{

class PointCloud
{
public:
  // Array of points
  std::vector<Point3D> Points;

  // Bounding box dimensions
  double XMin, YMin, XMax, YMax;

  // Create empty point cloud
  PointCloud() : XMin(0), YMin(0), XMax(0), YMax(0) {}
};

std::ostream &operator<<(std::ostream &stream, const PointCloud &pointCloud)
{
  stream << "Point cloud with " << pointCloud.Points.size()
         << " points on domain [" << pointCloud.XMin << ", " << pointCloud.XMax
         << "] x [" << pointCloud.YMin << ", " << pointCloud.YMax << "]";
  return stream;
}

} // namespace VirtualCity

#endif
