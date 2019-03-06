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


};

}

#endif
