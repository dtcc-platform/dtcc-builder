// Representation of a 2.5D building.
// Copyright (C) 2019-2020 Anders Logg, Anton J Olsson.
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

  /// ID (index) of building when parsed from shapefile
  int SHPFileID{};

  /// FNR of property containing building
  size_t PropertyFNR{};

  /// UUID of property containing building
  std::string PropertyUUID;

  /// Footprint of building (a 2D polygon)
  Polygon Footprint{};

  /// UUID of property containing building
  std::string BaseAreaID{};

  /// Height of building relative to ground height
  double Height{};

  /// Height of ground at site of building
  double GroundHeight{};

  /// Array of ground points (from point cloud) close to building footprint
  std::vector<Point3D> GroundPoints{};

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

  /// Check if building is empty
  bool Empty() const { return Footprint.Vertices.empty(); }

  /// Clear all data
  void Clear()
  {
    UUID = "";
    SHPFileID = 0;
    Footprint.Vertices.clear();
    Height = 0.0;
    GroundHeight = 0.0;
    GroundPoints.clear();
    RoofPoints.clear();
  }

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