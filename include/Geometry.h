// This class provides a collection of geometry utility functions.
// Anders Logg 2019

#ifndef VC_GEOMETRY_H
#define VC_GEOMETRY_H

#include <vector>
#include <cmath>
#include <limits>

#include "Point.h"
#include "Polygon.h"

namespace VirtualCity
{

class Geometry
{
public:

    // Compute dot product (2D)
    static double Dot2D(const Point2D& u, const Point2D& v)
    {
        return u.x * v.x + u.y * v.y;
    }

    // Compute dot product (3D)
    static double Dot3D(const Point3D& u, const Point3D& v)
    {
        return u.x * v.x + u.y * v.y + u.z * v.z;
    }

    // Compute distance between points (2D)
    static double Distance2D(const Point2D& p, const Point2D& q)
    {
        return std::sqrt(SquaredDistance2D(p, q));
    }

    // Compute distance between segment (p0, p1) and point q (2D)
    static double Distance2D(const Point2D& p0, const Point2D& p1,
                             const Point2D& q)
    {
        return std::sqrt(SquaredDistance2D(p0, p1, q));
    }

    // Compute distance between polygon and point (2D)
    static double Distance2D(const Polygon& polygon, const Point2D& p)
    {
        return std::sqrt(SquaredDistance2D(polygon, p));
    }

    // Compute distance between polygons (2D)
    static double Distance2D(const Polygon& polygon0, const Polygon& polygon1)
    {
        return std::sqrt(SquaredDistance2D(polygon0, polygon1));
    }

    // Compute distance between points (3D)
    static double Distance3D(const Point3D& p, const Point3D& q)
    {
        return std::sqrt(SquaredDistance3D(p, q));
    }

    // Compute squared distance between points (2D)
    static double SquaredDistance2D(const Point2D& p, const Point2D& q)
    {
        const double dx = p.x - q.x;
        const double dy = p.y - q.y;
        return dx * dx + dy * dy;
    }

    // Compute squared distance between segment (p0, p1) and point q (2D)
    static double SquaredDistance2D(const Point2D& p0, const Point2D& p1,
                                    const Point2D& q)
    {
        // Project point to line
        const Point2D u = q - p0;
        const Point2D v = p1 - p0;
        const Point2D p = p0 + v * (Dot2D(u, v) / v.SquaredMagnitude());

        // Check whether projected point is inside segment. Check either
        // x or y coordinates depending on which is largest (most stable)
        const bool inside =
            std::abs(v.x) > std::abs(v.y) ?
            std::min(p0.x, p1.x) <= p.x && p.x <= std::max(p0.x, p1.x) :
            std::min(p0.y, p1.y) <= p.y && p.y <= std::max(p0.y, p1.y);

        // Use distance to projection if inside
        if (inside)
            return (q - p).SquaredMagnitude();

        // Otherwise use distance to closest end point
        const double d0 = (q - p0).SquaredMagnitude();
        const double d1 = (q - p1).SquaredMagnitude();
        return std::min(d0, d1);
    }

    // Compute squared distance between polygon and point(2D)
    static double SquaredDistance2D(const Polygon& polygon, const Point2D& p)
    {
        // Check if point is contained in polygon
        if (PolygonContains2D(polygon, p))
            return 0.0;

        // If not, compute minimal squared distance to all segments
        double d2Min = std::numeric_limits<double>::max();
        for (size_t i = 0; i < polygon.Points.size(); i++)
        {
            Point2D p0 = polygon.Points[i];
            Point2D p1 = polygon.Points[(i + 1) % polygon.Points.size()];
            d2Min = std::min(d2Min, SquaredDistance2D(p0, p1, p));
        }

        return d2Min;
    }

    // Compute squared distance between polygons (2D)
    static double SquaredDistance2D(const Polygon& polygon0,
                                    const Polygon& polygon1)
    {
        double d2Min = std::numeric_limits<double>::max();

        // Check all vertices in first polygon
        for (auto const & p : polygon0.Points)
            d2Min = std::min(d2Min, SquaredDistance2D(polygon1, p));

        // Check all vertices in second polygon
        for (auto const & p : polygon1.Points)
            d2Min = std::min(d2Min, SquaredDistance2D(polygon0, p));

        return d2Min;
    }

    // Compute squared distance between points (3D)
    static double SquaredDistance3D(const Point3D& p, const Point3D& q)
    {
        const double dx = p.x - q.x;
        const double dy = p.y - q.y;
        const double dz = p.z - q.z;
        return dx * dx + dy * dy + dz * dz;
    }

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

    // Compute signed determinant of polygon (2D)
    static double PolygonDeterminant2D(const Polygon& polygon)
    {
        double sum = 0.0;
        for (size_t i = 0; i < polygon.Points.size(); i++)
        {
            Point2D p0 = polygon.Points[i];
            Point2D p1 = polygon.Points[(i + 1) % polygon.Points.size()];
            sum += (p1.x - p0.x) * (p1.y + p0.y);
        }
        return sum;
    }

    // Compute orientation of polygon (0 = counter-clockwise, 1 = clockwise)
    static size_t PolygonOrientation2D(const Polygon& polygon)
    {
        return PolygonDeterminant2D(polygon) < 0 ? 0 : 1;
    }

    // Compute area of polygon (2D)
    static double PolygonArea(const Polygon& polygon)
    {
        return 0.5 * std::abs(PolygonDeterminant2D(polygon));
    }

    // Compute center of polygon (2D)
    static Point2D PolygonCenter2D(const Polygon& polygon)
    {
        Point2D c;
        for (auto const & p : polygon.Points)
            c += p;
        c /= polygon.Points.size();
        return c;
    }

    // Compute radius of polygon relative to center
    static double PolygonRadius2D(const Polygon& polygon, const Point2D& center)
    {
        double r2max = 0.0;
        for (auto const & p : polygon.Points)
        {
            const double r2 = SquaredDistance2D(p, center);
            if (r2 > r2max)
                r2max = r2;
        }
        return std::sqrt(r2max);
    }

    // Check whether polygon contains point
    static bool PolygonContains2D(const Polygon& polygon, const Point2D& p)
    {
        // Compute total quadrant relative to polygon. If the point
        // is inside the polygon, the angle should be 4 (or -4).
        return Geometry::QuadrantAngle2D(p, polygon.Points) != 0;
    }

};

}

#endif
