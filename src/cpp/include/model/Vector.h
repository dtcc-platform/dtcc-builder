// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_VECTOR_H
#define DTCC_VECTOR_H

#include <cmath>

#include "Logging.h"
#include "Point.h"

namespace DTCC_BUILDER
{

/// Vector2D represents a Euclidean 2D vector
/// with the standard vector operations.
class Vector2D : public Printable
{
public:
  /// First component
  double x{};

  /// Second component
  double y{};

  /// Create zero vector
  Vector2D() = default;
  virtual ~Vector2D() {} // make the destructor virtual

  /// Create vector with given components.
  ///
  /// @param x First component
  /// @param y Second component
  Vector2D(double x, double y) : x(x), y(y) {}

  /// Create vector between origin and point (conversion from point).
  ///
  /// @param p The point
  explicit Vector2D(const Point2D &p) : x(p.x), y(p.y) {}

  /// Create vector between points.
  ///
  /// @param p First point
  /// @param q Second point
  explicit Vector2D(const Point2D &p, const Point2D &q)
      : x(q.x - p.x), y(q.y - p.y)
  {
  }

  /// Return point at origin + vector (conversion to point).
  ///
  /// @return Point at origin + vector.
  operator Point2D() const { return Point2D(x, y); }

  // FIXME: This class requires documentation

  Vector2D operator+(const Vector2D &p) const
  {
    Vector2D q(x + p.x, y + p.y);
    return q;
  }

  Vector2D operator+=(const Vector2D &p)
  {
    x += p.x;
    y += p.y;
    return *this;
  }

  Vector2D operator-(const Vector2D &p) const
  {
    return Vector2D{x - p.x, y - p.y};
  }

  Vector2D operator-() const { return Vector2D{-x, -y}; }

  Vector2D operator-=(const Vector2D &p)
  {
    x -= p.x;
    y -= p.y;
    return *this;
  }

  Vector2D operator*(double a) const
  {
    Vector2D q(x * a, y * a);
    return q;
  }

  Vector2D operator*=(double a)
  {
    x *= a;
    y *= a;
    return *this;
  }

  Vector2D operator/(double a) const
  {
    Vector2D q(x / a, y / a);
    return q;
  }

  Vector2D operator/=(double a)
  {
    x /= a;
    y /= a;
    return *this;
  }

  double Dot(const Vector2D &p) const { return x * p.x + y * p.y; }

  double AngleBetween(const Vector2D &p) const
  {
    return acos((Dot(p) / (Magnitude() * p.Magnitude())));
  }

  double Magnitude() const { return sqrt(SquaredMagnitude()); }

  double SquaredMagnitude() const { return x * x + y * y; }

  void Normalize() { (*this) /= Magnitude(); }

  /// Pretty-print
  std::string __str__() const override
  {
    return "(" + str(x) + ", " + str(y) + ")";
  }
};

/// Vector3D represents a Euclidean 3D vector
/// with the standard vector operations.
class Vector3D : public Printable
{
public:
  /// First component
  double x{};

  /// Second component
  double y{};

  /// Third component
  double z{};

  /// Create zero vector
  Vector3D() = default;
  virtual ~Vector3D() {} // make the destructor virtual
  /// Create vector with given components.
  ///
  /// @param x First component
  /// @param y Second component
  /// @param z Third component
  Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}

  /// Create vector between origin and point (conversion from point).
  ///
  /// @param p The point
  explicit Vector3D(const Point3D &p) : x(p.x), y(p.y), z(p.z) {}

  /// Create vector between points.
  ///
  /// @param p First point
  /// @param q Second point
  Vector3D(const Point3D &p, const Point3D &q)
      : x(q.x - p.x), y(q.y - p.y), z(q.z - p.z)
  {
  }

  /// Return point at origin + vector (conversion to point).
  ///
  /// @return Point at origin + vector.
  operator Point3D() const { return Point3D(x, y, z); }

  // FIXME: This class requires documentation

  Vector3D operator+(const Vector3D &p) const
  {
    Vector3D q(x + p.x, y + p.y, z + p.z);
    return q;
  }

  Vector3D operator+=(const Vector3D &p)
  {
    x += p.x;
    y += p.y;
    z += p.z;
    return *this;
  }

  Vector3D operator-(const Vector3D &p) const
  {
    return Vector3D{x - p.x, y - p.y, z - p.z};
  }

  Vector3D operator-() const { return Vector3D{-x, -y, -z}; }

  Vector3D operator-=(const Vector3D &p)
  {
    x -= p.x;
    y -= p.y;
    z -= p.z;
    return *this;
  }

  Vector3D operator*(double a) const
  {
    Vector3D q(x * a, y * a, z * a);
    return q;
  }

  Vector3D operator*=(double a)
  {
    x *= a;
    y *= a;
    z *= a;
    return *this;
  }

  Vector3D operator/(double a) const
  {
    Vector3D q(x / a, y / a, z / a);
    return q;
  }

  Vector3D operator/=(double a)
  {
    x /= a;
    y /= a;
    z /= a;
    return *this;
  }

  double Dot(const Vector3D &p) const { return x * p.x + y * p.y + z * p.z; }

  double Dot(const Point3D &p) const { return x * p.x + y * p.y + z * p.z; }

  Vector3D Cross(const Vector3D &p) const
  {
    return Vector3D{y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x};
  }

  double AngleBetween(const Vector3D &p) const
  {
    double a = Dot(p) / (Magnitude() * p.Magnitude());
    if (a > 1) // can happen due to rounding errors
      a = 1;
    return acos(a);
  }

  double Magnitude() const { return sqrt(SquaredMagnitude()); }

  double SquaredMagnitude() const { return x * x + y * y + z * z; }

  /// Pretty-print
  std::string __str__() const override
  {
    return "(" + str(x) + ", " + str(y) + ", " + str(z) + ")";
  }
};

// Note: We allow a minimal set of algebra for points such as
// translation by a vector. These need to be declared outside
// of the Point class to avoid circular includes.

/// Translate point by given vector.
///
/// @param p The point
/// @param v Translation vector
/// @return Translated point
Point2D operator+(const Point2D &p, const Vector2D &v)
{
  return {p.x + v.x, p.y + v.y};
}

/// Translate point by negative of given vector.
///
/// @param p The point
/// @param v Translation vector
/// @return Translated point
Point2D operator-(const Point2D &p, const Vector2D &v)
{
  return {p.x - v.x, p.y - v.y};
}

/// Translate point by given vector.
///
/// @param p The point
/// @param v Translation vector
/// @return Translated point
Point2D operator+=(Point2D &p, const Vector2D &v)
{
  p.x += v.x;
  p.y += v.y;
  return p;
}

/// Translate point by negative of given vector.
///
/// @param p The point
/// @param v Translation vector
/// @return Translated point
Point2D operator-=(Point2D &p, const Vector2D &v)
{
  p.x -= v.x;
  p.y -= v.y;
  return p;
}

/// Translate point by given vector.
///
/// @param p The point
/// @param v Translation vector
/// @return Translated point
Point3D operator+(const Point3D &p, const Vector3D &v)
{
  return {p.x + v.x, p.y + v.y, p.z + v.z};
}

/// Translate point by negative of given vector.
///
/// @param p The point
/// @param v Translation vector
/// @return Translated point
Point3D operator-(const Point3D &p, const Vector3D &v)
{
  return {p.x - v.x, p.y - v.y, p.z - v.z};
}

/// Translate point by given vector.
///
/// @param p The point
/// @param v Translation vector
/// @return Translated point
Point3D operator+=(Point3D &p, const Vector3D &v)
{
  p.x += v.x;
  p.y += v.y;
  p.z += v.z;
  return p;
}

/// Translate point by negative of given vector.
///
/// @param p The point
/// @param v Translation vector
/// @return Translated point
Point3D operator-=(Point3D &p, const Vector3D &v)
{
  p.x -= v.x;
  p.y -= v.y;
  p.z -= v.z;
  return p;
}

} // namespace DTCC_BUILDER

#endif
