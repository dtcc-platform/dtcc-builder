// Representation of a 2.5D city model.
// Copyright (C) 2019 Anders Logg.

#ifndef VC_CITY_MODEL_H
#define VC_CITY_MODEL_H

#include <string>
#include <iomanip>
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
            // Check whether point is inside building
            if (Geometry::PolygonContains2D(Buildings[i].Footprint, p))
                return i;
        }

        // Point not inside a building
        return -1;
    }

    // Compute center of city model
    Point2D Center() const
    {
        Point2D c;
        size_t numPoints = 0;
        for (auto const & building : Buildings)
        {
            for (auto const & p : building.Footprint.Points)
            {
                c += p;
                numPoints += 1;
            }
        }
        c /= numPoints;
        return c;
    }

    // Compute radius of city model (relative to center)
    double Radius(const Point2D& center) const
    {
        double r2max = 0.0;
        for (auto const & building : Buildings)
        {
            for (auto const & p : building.Footprint.Points)
            {
                const double r2 = Geometry::SquaredDistance2D(p, center);
                if (r2 > r2max)
                    r2max = r2;
            }
        }
        return std::sqrt(r2max);
    }

};

std::ostream& operator<<(std::ostream& stream, const CityModel& cityModel)
{
    const Point2D c = cityModel.Center();
    const double r = cityModel.Radius(c);
    stream << "CityModel with " << cityModel.Buildings.size()
           << " buildings and radius R = "
           << std::fixed << std::setprecision(2) << r;
}

}

#endif
