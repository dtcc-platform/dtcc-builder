// Point cloud (list of 3D points).
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

    // List of points
    std::vector<Point3D> Points;

};

}

#endif
