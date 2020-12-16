// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_BUILDING_H
#define DTCC_BUILDING_H

#include <string>
#include <vector>

#include "Logging.h"
#include "Polygon.h"
#include "Vector.h"

namespace DTCC
{

class Building : public Printable
{
public:
  // Building UUID
  std::string UUID;

  // Shapefile ID
  int SHPFileID{};

  // Uncomment for debugging
  // int debugID;

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

  /// Pretty-print
  std::string __str__() const override
  {
    return "Building with UUID " + UUID + " and height " + str(Height) +
           " and ground height " + str(GroundHeight);
  }
};

} // namespace DTCC

#endif
