// Height map (Digital Surface Model, DSM)
// Copyright (C) 2019 Anders Logg.

#ifndef VC_HEIGHT_MAP_H
#define VC_HEIGHT_MAP_H

#include <vector>

#include "Point.h"
namespace VirtualCity
{

class HeightMap
{
public:

    // Grid dimensions
    double XMin, YMin, XMax, YMax;

    // Number of grid points
    size_t XSize, YSize;

    // Grid data (flattened row-major starting at (XMin, YMin))
    std::vector<double> GridData;

    // Create empty height map
    HeightMap(double xMin, double yMin,
              double xMax, double yMax,
              double resolution)
        : XMin(xMin), YMin(yMin), XMax(xMax), YMax(yMax)
    {
        // Initialize grid data
        XSize = (XMax - XMin) / resolution + 1;
        YSize = (YMax - YMin) / resolution + 1;
        GridData.resize(XSize * YSize);
        std::fill(GridData.begin(), GridData.end(), 0.0);

        // Compute grid size
        hx = (XMax - XMin) / (XSize - 1);
        hy = (YMax - YMin) / (YSize - 1);
    }

    // Return height (z) at 2D point p
    double operator() (const Point2D& p) const
    {
        return (*this)(p.x, p.y);
    }

    // Return height (z) at 2D point (x, y)
    double operator() (double x, double y) const
    {
        return 0.0;
    }

    // Map index to coordinate
    Point2D Index2Coordinate(size_t i) const
    {
        const size_t ix = i % XSize;
        const size_t iy = i / YSize;
        return Point2D(XMin + ix * hx, YMin + iy * hy);
    }

    // Map coordinate to index (closest point)
    size_t Coordinate2Index(const Point2D& p) const
    {
        long int ix = std::lround((p.x - XMin) / hx);
        long int iy = std::lround((p.y - YMin) / hy);
        if (ix < 0) ix = 0;
        if (iy < 0) iy = 0;
        if (ix >= XSize) ix = XSize - 1;
        if (iy >= YSize) iy = YSize - 1;
        return iy * XSize + ix;
    }

private:

    // Grid size
    double hx, hy;

};

}

#endif
