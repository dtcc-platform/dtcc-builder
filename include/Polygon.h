// Polygon (array of 2D points).
// Copyright (C) 2019 Anders Logg.

#ifndef VC_POLYGON_H
#define VC_POLYGON_H

#include <vector>
#include <cmath>

#include "Point.h"
#include "Geometry.h"

namespace VirtualCity
{

class Polygon
{
public:

    // Array of points
    std::vector<Point2D> Points;

    // Create empty polygon
    Polygon() {}

    // Compute center of polygon
    Point2D Center() const
    {
        Point2D c;
        for (auto const & p : Points)
            c += p;
        c /= Points.size();
        return c;
    }

    // Compute radius of polygon
    double Radius(const Point2D& center) const
    {
        double r2max = 0.0;
        for (auto const & p : Points)
        {
            const double r2 = Geometry::SquaredDistance2D(p, center);
            if (r2 > r2max)
                r2max = r2;
        }
        return std::sqrt(r2max);
    }

    // Check whether polygon contains point
    bool Contains(const Point2D& p) const
    {
        // Compute total quadrant relative to polygon. If the point
        // is inside the polygon, the angle should be 4 (or -4).
        return Geometry::QuadrantAngle2D(p, Points) != 0;
    }

};

}

#endif
