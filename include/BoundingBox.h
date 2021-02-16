// Copyright (C) 2020-2021 Anders Logg, Anton J Olsson
// Licensed under the MIT License

#ifndef DTCC_BOUNDING_BOX_H
#define DTCC_BOUNDING_BOX_H

#include <cassert>
#include <limits>

#include "Point.h"
#include "Polygon.h"
#include "Logging.h"

namespace DTCC
{

  /// BoundingBox2D represents a 2D bounding box defined by two 2D points.
  class BoundingBox2D : public Printable
  {
  public:

    /// First ("lower left") corner
    Point2D P;

    /// Second ("upper right") corner
    Point2D Q;

    /// Create empty bounding box
    BoundingBox2D() = default;

    /// Create bounding box for given points.
    ///
    /// @param p First ("lower left") corner
    /// @param q Second ("upper right") corner
    BoundingBox2D(const Point2D& p, const Point2D& q) : P(p), Q(q)
    {
      assert(p.x <= q.x);
      assert(p.y <= q.y);
    }

    /// Create bounding box of vector of points.
    ///
    /// @param points Vector if points
    /// @param margin Margin to use for bounding box
    explicit BoundingBox2D(const std::vector<Point2D> &points,
                           double margin = 0.0)
    {
      constexpr double max = std::numeric_limits<double>::max();
      P.x = P.y = max;
      Q.x = Q.y = -max;
      for (const auto &p : points)
      {
        P.x = std::min(P.x, p.x);
        P.y = std::min(P.y, p.y);
        Q.x = std::max(Q.x, p.x);
        Q.y = std::max(Q.y, p.y);
      }
      P.x -= margin;
      P.y -= margin;
      Q.x += margin;
      Q.y += margin;
    }

    /// Return area of bounding box
    double Area() const { return (Q.x - P.x) * (Q.y - P.y); }

    /// Create bounding box of polygon.
    ///
    /// @param polygon Polygon
    explicit BoundingBox2D(const Polygon &polygon)
        : BoundingBox2D(polygon.Vertices)
    {
    }

    /// Pretty-print
    std::string __str__() const override
    {
      return
        "[" + str(P.x) + ", " + str(Q.x) + "] x " +
        "[" + str(P.y) + ", " + str(Q.y) + "]";
    }

  };

  /// BoundingBox3D represents a 3D bounding box defined by two 3D points.
  class BoundingBox3D : public Printable
  {
  public:

    /// First ("lower left") corner
    Point3D P{};

    /// Second ("upper right") corner
    Point3D Q{};

    /// Create empty bounding box
    BoundingBox3D() = default;

    /// Create bounding box defined by given points.
    ///
    /// @param p First ("lower left") corner
    /// @param q Second ("upper right") corner
    BoundingBox3D(const Point3D& p, const Point3D& q) : P(p), Q(q)
    {
      assert(p.x <= q.x);
      assert(p.y <= q.y);
      assert(p.z <= q.z);
    }

    /// Create bounding box of vector of points.
    ///
    /// @param points Vector if points
    /// @param margin Margin to use for bounding box
    explicit BoundingBox3D(const std::vector<Point3D> &points,
                           double margin = 0.0)
    {
      constexpr double max = std::numeric_limits<double>::max();
      P.x = P.y = P.z = max;
      Q.x = Q.y = Q.z = -max;
      for (const auto &p : points)
      {
        P.x = std::min(P.x, p.x);
        P.y = std::min(P.y, p.y);
        P.z = std::min(P.z, p.z);
        Q.x = std::max(Q.x, p.x);
        Q.y = std::max(Q.y, p.y);
        Q.z = std::max(Q.z, p.z);
      }
      P.x -= margin;
      P.y -= margin;
      P.z -= margin;
      Q.x += margin;
      Q.y += margin;
      Q.z += margin;
    }

    /// Return volume of bounding box
    double Volume() const { return (Q.x - P.x) * (Q.y - P.y) * (Q.z - P.z); }

    /// Pretty-print
    std::string __str__() const override
    {
      return
        "[" + str(P.x) + ", " + str(Q.x) + "] x " +
        "[" + str(P.y) + ", " + str(Q.y) + "] x " +
        "[" + str(P.z) + ", " + str(Q.z) + "]";
    }

  };

}

#endif
