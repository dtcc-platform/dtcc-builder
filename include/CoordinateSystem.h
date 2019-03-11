// Coordinate system transformation.
// Copyright (C) 2019 Anders Logg.

#ifndef VC_COORDINATE_SYSTEM_H
#define VC_COORDINATE_SYSTEM_H

#include <proj.h>
#include <sstream>

#include "Point.h"

namespace VirtualCity
{

class CoordinateSystem
{
public:

    // Example usage: Transform(p, "EPSG:4326", "EPSG:3006")

    // Transform point from between coordinate systems U and V
    static Point2D Transform(const Point2D& p, std::string U, std::string V)
    {
        // Create projection string
        std::stringstream s;
        s << "cs2cs +init=" << U << " +to +init=" << V;

        // Create projection
        PJ* P = proj_create (PJ_DEFAULT_CTX, s.str());

        // Transform coordinate
        PJ_COORD a = proj_coord (proj_torad(p.x), proj_torad(p.y), 0, 0);
        PJ_COORD b = proj_trans(P, PJ_FWD, a);

        // Clean up
        proj_destroy(P);
    }

};

}
