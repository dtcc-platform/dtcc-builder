// Representation of a 2.5D building.
// Copyright (C) 2019 Anders Logg.

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

  // Building UUID
  std::string UUID;
  int debugID;

  // Building height (above ground)
  double Height{};

  // Create empty building
  Building() {}
};

} // namespace DTCC

#endif
