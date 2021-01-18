// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_POINT_CLOUD_H
#define DTCC_POINT_CLOUD_H

#include <vector>

#include "BoundingBox.h"
#include "BoundingBoxTree.h"
#include "Color.h"
#include "Logging.h"
#include "Vector.h"

namespace DTCC
{

  class PointCloud : public Printable
  {
  public:

    /// Array of points
    std::vector<Vector3D> Points{};

    /// Array of colors (one per point)
    std::vector<Color> Colors{};

    /// Array of classifications (one per point)
    std::vector<uint8_t> Classification{};

    /// Bounding box
    BoundingBox2D BoundingBox{};

    /// Create empty point cloud
    PointCloud() = default;

    /// Return density of point cloud (points per square meter)
    double Density() const
    {
      return static_cast<double>(Points.size()) / BoundingBox.Area();
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

    /// Set a default color to all points
    void InitClassification(uint8_t c)
    {
      Classification.clear();
      for (size_t i = 0; i < Points.size(); i++)
      {
        Classification.push_back(c);
      }

    }

    /// Build search tree (bounding box tree), required for search queries.
    ///
    /// @param rebuild Force rebuild of existing tree if set
    void BuildSearchTree(bool rebuild = false) const
    {
      // Skip if already built or force rebuild
      if (!bbtree.Empty() && !rebuild)
      {
        Info("Search tree already built; set rebuild flag to force rebuild.");
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
      Progress(str(bbtree));
    }

    /// Clear all data
    void Clear()
    {
      Points.clear();
      Colors.clear();
      Classification.clear();
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
    // Bounding box tree used for search queries
    mutable BoundingBoxTree2D bbtree;
  };

} // namespace DTCC

#endif
