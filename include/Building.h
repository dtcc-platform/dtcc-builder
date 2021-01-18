// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_BUILDING_H
#define DTCC_BUILDING_H

#include <string>
#include <vector>

#include "Logging.h"
#include "Point.h"
#include "Polygon.h"
#include "Vector.h"

namespace DTCC
{

/// Building represents a building defined by a footprint and a height.
/// This means that a Building is currently an LOD1 representation.
class Building : public Printable
{
public:
  /// Universally unique identifier of building
  std::string UUID{};

  /// ID (index) of bulding when parsed from shapefile
  int SHPFileID{};

  /// Footprint of building (a 2D polygon)
  Polygon Footprint{};

  /// Height of building relative to ground height
  double Height{};

  /// Height of ground at site of building
  double GroundHeight{};

  /// Array of roof points (from point cloud) inside building footprint
  std::vector<Point3D> RoofPoints{};

  /// Create empty building
  Building() = default;

  // Uncomment for debugging
  // int debugID;

  /// Set new origin (subtract offset)
  void SetOrigin(const Point2D &origin) { Footprint.SetOrigin(origin); }

  /// Return minimum absolute height of building (equal to ground height).
  ///
  /// @return Minimum absolute height of building
  double MinHeight() const { return GroundHeight; }

  /// Return maximum absolute height of building
  ///
  /// @return Maximum absolute height of building
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
