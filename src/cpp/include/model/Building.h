// Representation of a 2.5D building.
// Copyright (C) 2019-2020 Anders Logg, Anton J Olsson.
// Licensed under the MIT License

#ifndef DTCC_BUILDING_H
#define DTCC_BUILDING_H

#include <string>
#include <vector>

#include "Logging.h"
#include "model/Point.h"
#include "model/Polygon.h"
#include "model/Vector.h"

namespace DTCC_BUILDER
{

// Error codes
enum BuildingError
{
  BUILDING_TOO_SMALL = 1 << 0,
  BUILDING_TOO_FEW_POINTS = 1 << 1,
  BUILDING_NO_ROOF_POINTS = 1 << 2,
  BUILDING_NO_GROUND_POINTS = 1 << 3,
  BUILDING_HEIGHT_TOO_LOW = 1 << 4,
  BUILDING_BAD_ASPECT_RATIO = 1 << 5,
  BUILDING_INSUFFICIENT_POINT_COVERAGE = 1 << 6
};

/// Building represents a building defined by a footprint and a height.
/// This means that a Building is currently an LOD1 representation.
class Building : public Printable
{
public:
  /// Universally unique identifier of building
  std::vector<std::string> AttachedUUIDs{};

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

  /// Error code for error found during while generating building
  size_t error = 0;

  /// Height of ground at site of building
  double GroundHeight{};

  /// Array of ground points (from point cloud) close to building footprint
  std::vector<Point3D> GroundPoints{};

  /// Array of roof points (from point cloud) inside building footprint
  std::vector<Point3D> RoofPoints{};

  /// Arrary of array of
  std::vector<std::vector<size_t>> RoofSegments{};

  /// Create empty building
  Building() = default;
  virtual ~Building() {} // make the destructor virtual

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

  /// Check if building is valid (at least 3 vertices)
  bool Valid() const { return Footprint.Vertices.size() >= 3; }

  /// Clear all data
  void Clear()
  {
    UUID = "";
    SHPFileID = 0;
    Footprint.Vertices.clear();
    Height = 0.0;
    GroundHeight = 0.0;
    error = 0;
    GroundPoints.clear();
    RoofPoints.clear();
    RoofSegments.clear();
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

} // namespace DTCC_BUILDER

#endif
