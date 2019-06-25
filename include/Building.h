// Representation of a 2.5D building.
// Copyright (C) 2019 Anders Logg.

#ifndef VC_BUILDING_H
#define VC_BUILDING_H

#include <vector>

#include "Point.h"
#include "Polygon.h"

namespace VirtualCity
{

class Building
{
public:
  // Building footprint (polygon)
  Polygon Footprint;

  // Building height (above ground)
  double Height;

  // Create empty building
  Building() : Height(0) {}
};

} // namespace VirtualCity

#endif
