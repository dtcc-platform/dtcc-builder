// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_POINT_CLOUD_H
#define DTCC_POINT_CLOUD_H

#include <set>
#include <vector>

#include "BoundingBox.h"
#include "BoundingBoxTree.h"
#include "Color.h"
#include "Logging.h"
#include "Point.h"
#include "Vector.h"

namespace DTCC_BUILDER
{

class PointCloud : public Printable
{
public:
  /// Array of points
  std::vector<Point3D> Points{};

  /// Array of Normals
  std::vector<Vector3D> Normals{};

  /// Array of colors (one per point)
  std::vector<Color> Colors{};

  /// Array of classifications (one per point)
  std::vector<uint8_t> Classifications{};

  /// Array of intensities (one per point)
  std::vector<uint16_t> Intensities{};

  /// Array of scan data (one per point)
  std::vector<uint8_t> ScanFlags{};

  /// Bounding box
  BoundingBox2D BoundingBox{};

  Point2D Origin = Point2D(0, 0);

  /// Check if point cloud is empty
  bool Empty() const { return Points.empty(); }

  /// Create empty point cloud
  PointCloud() = default;

  /// Return density of point cloud (points per square meter)
  double Density() const
  {
    return static_cast<double>(Points.size()) / BoundingBox.Area();
  }

  void CalculateBoundingBox() { BoundingBox = BoundingBox2D(Points); }

  /// Set new origin (subtract offset). Note that this only
  /// affects the x and y coordinates (z unaffected).
  void SetOrigin(const Point2D &origin)
  {
    info("PointCloud: Setting new origin to " + str(origin));
    const Vector2D v2D{origin.x, origin.y};
    const Vector3D v3D{origin.x, origin.y, 0.0};
    for (auto &p : Points)
      p -= v3D;
    BoundingBox.P -= v2D;
    BoundingBox.Q -= v2D;
    Origin.x = origin.x;
    Origin.y = origin.y;
  }

  /// set a default color to all points
  void InitColors(const Color &c)
  {
    Colors.clear();
    for (size_t i = 0; i < Points.size(); i++)
    {
      Colors.push_back(c);
    }
  }

  /// Set a default classification to all points
  void InitClassifications(uint8_t c)
  {
    Classifications.clear();
    for (size_t i = 0; i < Points.size(); i++)
      Classifications.push_back(c);
  }

  /// Build search tree (bounding box tree), required for search queries.
  ///
  /// @param rebuild Force rebuild of existing tree if set
  void BuildSearchTree(bool rebuild = false) const
  {
    // Skip if already built or force rebuild
    if (!bbtree.Empty() && !rebuild)
    {
      info("Search tree already built; set rebuild flag to force rebuild.");
      return;
    }

    // Create 2D bounding boxes for all points
    std::vector<BoundingBox2D> bboxes;
    for (const auto &p3D : Points)
    {
      const Point2D p2D(p3D.x, p3D.y);
      BoundingBox2D bbox(p2D, p2D);
      bboxes.push_back(bbox);
    }

    // Build bounding box tree
    bbtree.Build(bboxes);
    debug(str(bbtree));
  }

  void BuildHasClassifications()
  {
    for (auto c : Classifications)
    {
      UsedClassifications.insert(c);
    }
  }

  bool HasClassification(uint8_t c) const

  {
    if (UsedClassifications.empty())
    {
      warning("Classification set is not built.");
      return false;
    }
    return UsedClassifications.find(c) != UsedClassifications.end();
  }

  /// Clear all data
  void Clear()
  {
    Points.clear();
    Normals.clear();
    Colors.clear();
    Classifications.clear();
    UsedClassifications.clear();
    Intensities.clear();
    ScanFlags.clear();
    bbtree.Clear();
  }

  /// Pretty-print
  std::string __str__() const override
  {
    return "Point cloud on " + str(BoundingBox) + " with " +
           str(Points.size()) + " points and density " + str(Density()) +
           " m^-2";
  }

private:
  /// Set of used point classifications
  std::set<uint8_t> UsedClassifications{};

  friend class CityModelGenerator;

  // Bounding box tree used for search queries
  mutable BoundingBoxTree2D bbtree;
};

} // namespace DTCC_BUILDER

#endif
