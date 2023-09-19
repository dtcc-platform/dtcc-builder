// Copyright (C) 2020-2021 Anders Logg, Anton J Olsson, Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_BOUNDING_BOX_H
#define DTCC_BOUNDING_BOX_H

#include <cassert>
#include <limits>

#include "Logging.h"
#include "model/Polygon.h"
#include "model/Vector.h"

namespace DTCC_BUILDER
{

/// BoundingBox2D represents a 2D bounding box defined by two 2D points.
class BoundingBox2D : public Printable
{
public:
  /// First ("lower left") corner
  Vector2D P;

  /// Second ("upper right") corner
  Vector2D Q;

  /// Create empty bounding box
  BoundingBox2D() = default;
  virtual ~BoundingBox2D() {} // make the destructor virtual

  /// Create bounding box for given points.
  ///
  /// @param p First ("lower left") corner
  /// @param q Second ("upper right") corner
  BoundingBox2D(const Vector2D &p, const Vector2D &q) : P(p), Q(q)
  {
    assert(p.x <= q.x);
    assert(p.y <= q.y);
  }

  /// Create bounding box of vector of points.
  ///
  /// @param points Vector if points
  /// @param margin Margin to use for bounding box
  explicit BoundingBox2D(const std::vector<Vector2D> &points,
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

  /// Create 2D bounding box of vector of 3D points.
  ///
  /// @param points Vector if points
  /// @param margin Margin to use for bounding box
  explicit BoundingBox2D(const std::vector<Vector3D> &points,
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

  /// Create bounding box of vector of polygons.
  ///
  /// @param points Vector of polygons
  /// @param margin Margin to use for bounding box
  explicit BoundingBox2D(const std::vector<Polygon> &polygons,
                         double margin = 0.0)
  {
    auto init_bb = BoundingBox2D(polygons.front());
    P = init_bb.P;
    Q = init_bb.Q;
    for (auto p : polygons)
    {
      this->union_with(BoundingBox2D(p));
    }
    P.x -= margin;
    P.y -= margin;
    Q.x += margin;
    Q.y += margin;
  }

  /// expand a bounding box to the union of it
  /// and a second bounding box
  void union_with(BoundingBox2D const &other)
  {
    P.x = std::min(P.x, other.P.x);
    P.y = std::min(P.y, other.P.y);
    Q.x = std::max(Q.x, other.Q.x);
    Q.y = std::max(Q.y, other.Q.y);
  }

  /// intersect bounding box with a second bounding box
  void intersect(BoundingBox2D const &other)
  {
    if (P.x >= other.Q.x || Q.x <= other.P.x || P.y >= other.Q.y ||
        Q.y <= other.P.y)
    {
      // empty BB
      P = Vector2D();
      Q = Vector2D();
    }
    else
    {
      P.x = std::max(P.x, other.P.x);
      Q.x = std::min(Q.x, other.Q.x);
      P.y = std::max(P.y, other.P.y);
      Q.y = std::min(Q.y, other.Q.y);
    }
  }

  /// Return area of bounding box
  double area() const { return (Q.x - P.x) * (Q.y - P.y); }

  /// Create bounding box of polygon.
  ///
  /// @param polygon Polygon
  explicit BoundingBox2D(const Polygon &polygon)
      : BoundingBox2D(polygon.vertices)
  {
  }

  /// Pretty-print
  std::string __str__() const override
  {
    return "[" + str(P.x) + ", " + str(Q.x) + "] x " + "[" + str(P.y) + ", " +
           str(Q.y) + "]";
  }
};

/// BoundingBox3D represents a 3D bounding box defined by two 3D points.
class BoundingBox3D : public Printable
{
public:
  /// First ("lower left") corner
  Vector3D P{};

  /// Second ("upper right") corner
  Vector3D Q{};

  /// Create empty bounding box
  BoundingBox3D() = default;
  virtual ~BoundingBox3D() {} // make the destructor virtual

  /// Create bounding box defined by given points.
  ///
  /// @param p First ("lower left") corner
  /// @param q Second ("upper right") corner
  BoundingBox3D(const Vector3D &p, const Vector3D &q) : P(p), Q(q)
  {
    assert(p.x <= q.x);
    assert(p.y <= q.y);
    assert(p.z <= q.z);
  }

  /// Create bounding box of vector of points.
  ///
  /// @param points Vector if points
  /// @param margin Margin to use for bounding box
  explicit BoundingBox3D(const std::vector<Vector3D> &points,
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

  /// expand a bounding box to the union of it
  /// and a second bounding box
  void union_with(BoundingBox3D other)
  {
    P.x = std::min(P.x, other.P.x);
    P.y = std::min(P.y, other.P.y);
    P.z = std::min(P.z, other.P.z);
    Q.x = std::max(Q.x, other.Q.x);
    Q.y = std::max(Q.y, other.Q.y);
    Q.z = std::max(Q.z, other.Q.z);
  }

  /// Return volume of bounding box
  double volume() const { return (Q.x - P.x) * (Q.y - P.y) * (Q.z - P.z); }

  /// Pretty-print
  std::string __str__() const override
  {
    return "[" + str(P.x) + ", " + str(Q.x) + "] x " + "[" + str(P.y) + ", " +
           str(Q.y) + "] x " + "[" + str(P.z) + ", " + str(Q.z) + "]";
  }
};

} // namespace DTCC_BUILDER

#endif
