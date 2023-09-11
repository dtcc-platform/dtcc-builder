// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_POINT_H
#define DTCC_POINT_H

#include "Logging.h"

namespace DTCC_BUILDER
{

/// Point2D represents a point in 2D Euclidean space.
class Point2D : public Printable
{
public:
  /// First coordinate
  double x{};

  /// Second coordinate
  double y{};

  /// Create point at origin
  Point2D() = default;
  virtual ~Point2D() {} // make the destructor virtual

  /// Create point with given coordinates.
  ///
  /// @param x First coordinate
  /// @param y Second coordinate
  Point2D(double x, double y) : x(x), y(y) {}

  double operator[](int idx) const
  {
    if (idx == 0)
      return x;
    if (idx == 1)
      return y;
    throw std::out_of_range("Out of bounds");
  }

  size_t size() const { return 2; }

  /// Pretty-print
  std::string __str__() const override
  {
    return "(" + str(x) + ", " + str(y) + ")";
  }
};

/// Point3D represents a point in 3D Euclidean space.
class Point3D : public Printable
{
public:
  /// First coordinate
  double x{};

  /// Second coordinate
  double y{};

  /// Third coordinate
  double z{};

  /// Create point at origin
  Point3D() = default;
  virtual ~Point3D() {}

  /// Create point with given coordinates.
  ///
  /// @param x First coordinate
  /// @param y Second coordinate
  /// @param z Third coordinate
  Point3D(double x, double y, double z) : x(x), y(y), z(z) {}

  double operator[](int idx) const
  {
    if (idx == 0)
      return x;
    if (idx == 1)
      return y;
    if (idx == 2)
      return z;
    throw std::out_of_range("Out of bounds");
  }

  size_t size() const { return 3; }

  /// Pretty-print
  std::string __str__() const override
  {
    return "(" + str(x) + ", " + str(y) + ", " + str(z) + ")";
  }
};

} // namespace DTCC_BUILDER

#endif
