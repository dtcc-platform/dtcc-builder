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
  // Building UUID
  std::string UUID;
  int debugID;
  int entityID;

  // Building footprint (polygon)
  Polygon Footprint{};

  // Building height (relative to ground)
  double Height{};

  // Ground height (absolute)
  double GroundHeight{};

  // Create empty building
  Building() = default;

  // Return minimum absolute height of building
  double MinHeight() const { return GroundHeight; }

  // Return maximum absolute height of building
  double MaxHeight() const { return GroundHeight + Height; }
};

} // namespace DTCC

#endif
