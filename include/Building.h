// Representation of a 2.5D building.
// Copyright (C) 2019 Anders Logg.
// Licensed under the MIT License

#ifndef DTCC_BUILDING_H
#define DTCC_BUILDING_H

#include <vector>

#include "Vector.h"
#include "Polygon.h"

namespace DTCC
{

class Building
{
public:
  // Building footprint (polygon)
  Polygon Footprint{};

  // Building height (relative to ground)
  double Height{};

  // Ground height (absolute)
  double GroundHeight{};

  // Create empty building
  Building() = default;
};

} // namespace DTCC

#endif
