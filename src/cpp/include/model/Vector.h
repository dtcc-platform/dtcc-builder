// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_VECTOR_H
#define DTCC_VECTOR_H

#include <cmath>

#include "Logging.h"
//#include "Point.h"

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

  /// Create vector between vectors.
  ///
  /// @param p First point
  /// @param q Second point
  explicit Vector2D(const Vector2D &p, const Vector2D &q)
      : x(q.x - p.x), y(q.y - p.y)
  {
  }

  // Operator merged from Point2D class
  double operator[](int idx) const
  {
    if (idx == 0)
      return x;
    if (idx == 1)
      return y;
    throw std::out_of_range("Out of bounds");
  }

  size_t size() const { return 2; }

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

  double dot(const Vector2D &p) const { return x * p.x + y * p.y; }

  double angle_between(const Vector2D &p) const
  {
    return acos((dot(p) / (magnitude() * p.magnitude())));
  }

  double magnitude() const { return sqrt(squared_magnitude()); }

  double squared_magnitude() const { return x * x + y * y; }

  void normalize() { (*this) /= magnitude(); }

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

  /// Create vector between points/vectors.
  ///
  /// @param p First point
  /// @param q Second point
  Vector3D(const Vector3D &p, const Vector3D &q)
      : x(q.x - p.x), y(q.y - p.y), z(q.z - p.z)
  {
  }

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

  double dot(const Vector3D &p) const { return x * p.x + y * p.y + z * p.z; }

  Vector3D cross(const Vector3D &p) const
  {
    return Vector3D{y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x};
  }

  double angle_between(const Vector3D &p) const
  {
    double a = dot(p) / (magnitude() * p.magnitude());
    if (a > 1) // can happen due to rounding errors
      a = 1;
    return acos(a);
  }

  double magnitude() const { return sqrt(squared_magnitude()); }

  double squared_magnitude() const { return x * x + y * y + z * z; }

  /// Pretty-print
  std::string __str__() const override
  {
    return "(" + str(x) + ", " + str(y) + ", " + str(z) + ")";
  }
};

} // namespace DTCC_BUILDER

#endif
