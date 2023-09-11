// Representation of a 2.5D building.
// Copyright (C) 2019-2020 Anders Logg, Anton J Olsson.
// Licensed under the MIT License

#ifndef DTCC_BUILDING_H
#define DTCC_BUILDING_H

#include <string>
#include <vector>

#include "Logging.h"
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
  std::vector<std::string> attached_uuids{};

  /// Universally unique identifier of building
  std::string uuid{};

  /// ID (index) of building when parsed from shapefile
  int shpfile_id{};

  /// FNR of property containing building
  size_t property_fnr{};

  /// uuid of property containing building
  std::string property_uuid;

  /// footprint of building (a 2D polygon)
  Polygon footprint{};

  /// uuid of property containing building
  std::string base_area_id{};

  /// height of building relative to ground height
  double height{};

  /// Error code for error found during while generating building
  size_t error = 0;

  /// height of ground at site of building
  double ground_height{};

  /// Array of ground points (from point cloud) close to building footprint
  std::vector<Vector3D> ground_points{};

  /// Array of roof points (from point cloud) inside building footprint
  std::vector<Vector3D> roof_points{};

  /// Arrary of array of
  std::vector<std::vector<size_t>> roof_segments{};

  /// Create empty building
  Building() = default;
  virtual ~Building() {} // make the destructor virtual

  // Uncomment for debugging
  // int debugID;

  /// Set new origin (subtract offset)
  void set_origin(const Vector2D &origin) { footprint.set_origin(origin); }

  /// Return minimum absolute height of building (equal to ground height).
  ///
  /// @return Minimum absolute height of building
  double min_height() const { return ground_height; }

  /// Return maximum absolute height of building
  ///
  /// @return Maximum absolute height of building
  double max_height() const { return ground_height + height; }

  /// Check if building is empty
  bool empty() const { return footprint.vertices.empty(); }

  /// Check if building is valid (at least 3 vertices)
  bool valid() const { return footprint.vertices.size() >= 3; }

  /// clear all data
  void clear()
  {
    uuid = "";
    shpfile_id = 0;
    footprint.vertices.clear();
    height = 0.0;
    ground_height = 0.0;
    error = 0;
    ground_points.clear();
    roof_points.clear();
    roof_segments.clear();
  }

  /// Pretty-print Building.
  /// \return Pretty-print string
  std::string __str__() const override
  {
    return "Building\n"
           "uuid: " +
           uuid + "\nProperty FNR: " + str(property_fnr) +
           "\nProperty uuid: " + property_uuid +
           "\nBase area ID: " + base_area_id +
           "\nMin height: " + std::to_string(min_height()) +
           "\nMax height: " + std::to_string(max_height()) +
           "\nFootprint: " + footprint.__str__();
  }
};

} // namespace DTCC_BUILDER

#endif
