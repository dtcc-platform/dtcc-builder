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
#include "Vector.h"

namespace DTCC_BUILDER
{

class PointCloud : public Printable
{
public:
  /// Array of points
  std::vector<Vector3D> points{};

  /// Array of normals
  std::vector<Vector3D> normals{};

  /// Array of colors (one per point)
  std::vector<Color> colors{};

  /// Array of classifications (one per point)
  std::vector<uint8_t> classifications{};

  /// Array of intensities (one per point)
  std::vector<uint16_t> intensities{};

  /// Array of scan data (one per point)
  std::vector<uint8_t> scan_flags{};

  /// Bounding box
  BoundingBox2D bounding_box{};

  Vector2D origin = Vector2D(0, 0);

  /// Check if point cloud is empty
  bool empty() const { return points.empty(); }

  /// Create empty point cloud
  PointCloud() = default;
  virtual ~PointCloud() {} // make the destructor virtual

  /// Return density of point cloud (points per square meter)
  double density() const
  {
    return static_cast<double>(points.size()) / bounding_box.area();
  }

  void calculate_bounding_box() { bounding_box = BoundingBox2D(points); }

  /// Set new origin (subtract offset). Note that this only
  /// affects the x and y coordinates (z unaffected).
  void set_origin(const Vector2D &origin)
  {
    info("PointCloud: Setting new origin to " + str(origin));
    const Vector2D v_2d{origin.x, origin.y};
    const Vector3D v_3d{origin.x, origin.y, 0.0};
    for (auto &p : points)
      p -= v_3d;
    bounding_box.P -= v_2d;
    bounding_box.Q -= v_2d;
    this->origin.x = origin.x;
    this->origin.y = origin.y;
  }

  /// set a default color to all points
  void init_colors(const Color &c)
  {
    colors.clear();
    for (size_t i = 0; i < points.size(); i++)
    {
      colors.push_back(c);
    }
  }

  /// Set a default classification to all points
  void init_classifications(uint8_t c)
  {
    classifications.clear();
    for (size_t i = 0; i < points.size(); i++)
      classifications.push_back(c);
  }

  /// build search tree (bounding box tree), required for search queries.
  ///
  /// @param rebuild Force rebuild of existing tree if set
  void build_search_tree(bool rebuild = false) const
  {
    // Skip if already built or force rebuild
    if (!bbtree.empty() && !rebuild)
    {
      info("Search tree already built; set rebuild flag to force rebuild.");
      return;
    }

    // Create 2D bounding boxes for all points
    std::vector<BoundingBox2D> bboxes;
    for (const auto &p_3d : points)
    {
      const Vector2D p_2d(p_3d.x, p_3d.y);
      BoundingBox2D bbox(p_2d, p_2d);
      bboxes.push_back(bbox);
    }

    // build bounding box tree
    bbtree.build(bboxes);
    debug(str(bbtree));
  }

  void build_has_classifications()
  {
    for (auto c : classifications)
    {
      used_classifications.insert(c);
    }
  }

  bool has_classification(uint8_t c) const

  {
    if (used_classifications.empty())
    {
      warning("Classification set is not built.");
      return false;
    }
    return used_classifications.find(c) != used_classifications.end();
  }

  /// clear all data
  void clear()
  {
    points.clear();
    normals.clear();
    colors.clear();
    classifications.clear();
    used_classifications.clear();
    intensities.clear();
    scan_flags.clear();
    bbtree.clear();
  }

  /// Pretty-print
  std::string __str__() const override
  {
    return "Point cloud on " + str(bounding_box) + " with " +
           str(points.size()) + " points and density " + str(density()) +
           " m^-2";
  }

private:
  /// Set of used point classifications
  std::set<uint8_t> used_classifications{};

  friend class CityModelGenerator;

  // Bounding box tree used for search queries
  mutable BoundingBoxTree2D bbtree;
};

} // namespace DTCC_BUILDER

#endif
