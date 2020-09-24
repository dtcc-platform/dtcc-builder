// Coordinate system transformation.
// Copyright (C) 2019 Anders Logg.

#ifndef DTCC_COORDINATE_SYSTEM_H
#define DTCC_COORDINATE_SYSTEM_H

#include <proj_api.h>
#include <sstream>

#include "Vector.h"

namespace DTCC
{

class CoordinateSystem
{
public:
  // Example usage: q = Transform(p, "epsg:4326", "epsg:3006")
  // User can use use https://mygeodata.cloud/cs2cs/ for validation

  // Developer note: We are currently using PROJ version 4 instead
  // of the newer version 5 since only version 4 is available in
  // Ubuntu 18.04 (LTS). The version 4 API is deprecated and we will
  // need to migrate in the future. See notes on migration:
  // https://proj4.org/development/migration.html

  // Transform point from between coordinate systems U and V
  static Vector2D Transform(const Vector2D &p, std::string U, std::string V)
  {
    // Set up projections
    const std::string s;
    const std::string sU = s + "+init=" + U;
    const std::string sV = s + "+init=" + V;
    const projPJ pjU = pj_init_plus(sU.c_str());
    const projPJ pjV = pj_init_plus(sV.c_str());

    // Check errors
    if (!pjU)
      throw std::runtime_error(s + "Unknown input projection: " + sU);
    if (!pjV)
      throw std::runtime_error(s + "Unknown output projection: " + sV);

    // Set input coordinates
    double x = p.x;
    double y = p.y;

    // Compute transformation
    pj_transform(pjU, pjV, 1, 1, &x, &y, NULL);

    // Set output coordinates
    Vector2D q(x, y);

    return q;
  }
};

} // namespace DTCC

#endif
