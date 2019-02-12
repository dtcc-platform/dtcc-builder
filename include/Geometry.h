// This class provides a collection of geometry utility functions.
// Anders Logg 2019

#ifndef VC_GEOMETRY_H
#define VC_GEOMETRY_H

#include <vector>

#include "Point.h"

namespace VirtualCity
{

class Geometry
{
public:

    // Compute quadrant angle of point p relative to polygon (2D)
    static int QuadrantAngle2D(const Point2D& p,
                               const std::vector<Point2D>& polygon)
    {
        // Compute angle to first vertex
        Point2D q0 = polygon[0];
        int v0 = QuadrantAngle2D(q0, p);

        // Sum up total angle
        int totalAngle = 0;
        for (int i = 1; i < polygon.size() + 1; i++)
        {
            // Compute angle increment
            Point2D q1 = polygon[i % polygon.size()];
            int v1 = QuadrantAngle2D(q1, p);
            int dv = v1 - v0;

            // Adjust angle increment for wrap-around
            if (dv == 3)
                dv = -1;
            else if (dv == -3)
                dv = 1;
            else if (dv == 2 || dv == -2)
            {
                double xx = q1.x - ((q1.y - p.y) * ((q0.x - q1.x) / (q0.y - q1.y)));
                if (xx > p.x)
                    dv = -dv;
            }

            // Add to total angle and update
            totalAngle += dv;
            q0 = q1;
            v0 = v1;
        }

        return totalAngle;
    }

    // Compute quadrant angle of point p relative to point q (2D)
    static int QuadrantAngle2D(const Point2D& p, const Point2D& q)
    {
        return ((p.x > q.x) ? ((p.y > q.y) ? 0 : 3) : ((p.y > q.y) ? 1 : 2));
    }

};

}

#endif
