// Representation of a 2.5D building.
// Copyright (C) 2019-2020 Anders Logg, Anton J Olsson.
// Licensed under the MIT License

#ifndef DTCC_BUILDING_H
#define DTCC_BUILDING_H

#include <vector>

#include "Vector.h"
#include "Polygon.h"

namespace DTCC
{

class Building : public Printable
{
public:
  /// Building's UUID
  std::string UUID;

  // Uncomment for debugging
  // int debugID;

  /// Shapefile ID
  int SHPFileID{};

  /// FNR of property containing building
  size_t PropertyFNR{};

  /// UUID of property containing building
  std::string PropertyUUID;

  /// Building footprint (polygon)
  Polygon Footprint{};

  /// UUID of property containing building
  std::string BaseAreaID{};

  /// Building height (relative to ground)
  double Height{};

  /// Ground height (absolute)
  double GroundHeight{};

  /// Create empty building.
  Building() = default;

  /// Return minimum absolute height of building.
  double MinHeight() const { return GroundHeight; }

  /// Return maximum absolute height of building.
  double MaxHeight() const { return GroundHeight + Height; }

  /// Pretty-print Building.
  /// \return Pretty-print string
  std::string __str__() const override
  {
    return "Building\n"
           "UUID: " +
           UUID + "\nProperty FNR: " + str(PropertyFNR) +
           "\nProperty UUID: " + PropertyUUID +
           "\nBase area ID: " + BaseAreaID +
           "\nMin height: " + std::to_string(MinHeight()) +
           "\nMax height: " + std::to_string(MaxHeight()) +
           "\nFootprint: " + Footprint.__str__();
  }
};

} // namespace DTCC

#endif
