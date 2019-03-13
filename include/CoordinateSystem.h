// Coordinate system transformation.
// Copyright (C) 2019 Anders Logg.

#ifndef VC_COORDINATE_SYSTEM_H
#define VC_COORDINATE_SYSTEM_H

#include <proj_api.h>
#include <sstream>

#include "Point.h"

namespace VirtualCity
{

class CoordinateSystem
{
public:

    // Example usage: Transform(p, "EPSG:4326", "EPSG:3006")
    // User can use use https://mygeodata.cloud/cs2cs/ for validation

    // Developer note: We are currently using PROJ version 4 instead
    // of the newer version 5 since only version 4 is available in
    // Ubuntu 18.04 (LTS). The version 4 API is deprecated and we will
    // need to migrate in the future. See notes on migration:
    // https://proj4.org/development/migration.html

    // Transform point from between coordinate systems U and V
    static Point2D Transform(const Point2D& p, std::string U, std::string V)
    {

        projPJ pj_merc = pj_init_plus("+proj=merc +ellps=clrk66 +lat_ts=33");
        if (!pj_merc)
            throw std::runtime_error("Error 1");
        projPJ pj_latlong = pj_init_plus("+proj=latlong +ellps=clrk66");
        if (!pj_latlong)
            throw std::runtime_error("Error 2");


        double x = p.x;
        double y = p.y;
        x *= DEG_TO_RAD;
        y *= DEG_TO_RAD;
        pj_transform(pj_latlong, pj_merc, 1, 1, &x, &y, NULL );
        printf("%.2f\t%.2f\n", x, y);

        return Point2D();

        // // Create projection string
        // std::stringstream s;
        // s << "cs2cs +init=" << U << " +to +init=" << V;

        // // Create projection
        // PJ* P = proj_create(PJ_DEFAULT_CTX, s.str());

        // // Transform coordinate
        // PJ_COORD a = proj_coord(proj_torad(p.x), proj_torad(p.y), 0, 0);
        // PJ_COORD b = proj_trans(P, PJ_FWD, a);

        // // Clean up
        // proj_destroy(P);
    }

};

}

#endif
