// Representation of a 2.5D city model.
// Copyright (C) 2019 Anders Logg.

#ifndef VC_CITY_MODEL_H
#define VC_CITY_MODEL_H

#include <vector>

#include "Point.h"
#include "Building.h"
#include "Geometry.h"

namespace VirtualCity
{

class CityModel
{
public:

    // List of buildings
    std::vector<Building> Buildings;

    // Find building containing point (inside footprint), returning -1
    // if the point is not inside any building.
    int FindBuilding(const Point2D& p) const
    {
        // Iterate over buildings
        for (size_t i = 0; i < Buildings.size(); i++)
        {
            // Compute total quadrant relative to polygon. If the point
            // is inside the polygon, the angle should be 4 (or -4).
            const int v = Geometry::QuadrantAngle2D(p, Buildings[i].Footprint);

            // Check if point is inside the domain
            if (v != 0)
                return i;
        }

        // Point not inside a building
        return -1;
    }

    // Compute domain center
    Point2D DomainCenter() const
    {
        Point2D c;
        double r;
        ComputeDomainSize(c, r);
        return c;
    }

    // Compute domain radius
    double DomainRadius() const
    {
        Point2D c;
        double r;
        ComputeDomainSize(c, r);
        return r;
    }

    // Compute domain size (center and radius)
    void ComputeDomainSize(Point2D& c, double& r) const
    {
        // The center and radius is computed by brute force and is O(N^2).
        // It can be more efficiently computed from the convex hull but
        // this step is likely relatively inexpensive.

        // Compute diameter and keep points at largest distance
        Point2D q0, q1;
        double d2max = 0.0;
        for (auto const & b0 : Buildings)
        {
            for (auto const & b1 : Buildings)
            {
                for (auto const & p0 : b0.Footprint)
                {
                    for (auto const & p1 : b1.Footprint)
                    {
                        const double d2 = Geometry::SquaredDistance2D(p0, p1);
                        if (d2 > d2max)
                        {
                            q0 = p0;
                            q1 = p1;
                            d2max = d2;
                        }
                    }
                }
            }
        }

        // Compute center and radius
        c = (q0 + q1) * 0.5;
        r = 0.5 * std::sqrt(d2max);
    }

};

std::ostream& operator<<(std::ostream& stream, const CityModel& cityModel)
{
    Point2D c;
    double r;
    cityModel.ComputeDomainSize(c, r);
    stream << "CityModel with " << cityModel.Buildings.size()
           << " buildings, radius R = " << r
           << " and center C = " << c << std::endl;
}

}

#endif
