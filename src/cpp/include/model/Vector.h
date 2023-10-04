// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_VECTOR_H
#define DTCC_VECTOR_H

#include <cmath>
#include <functional>

#include "Logging.h"
// #include "Point.h"

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

  bool operator==(const Vector3D &p) const
  {
    return (x == p.x && y == p.y && z == p.z);
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

  Vector3D rotate(const Vector3D &axis, double angle) const
  {
    double c = cos(angle);
    double s = sin(angle);
    double t = 1 - c;
    double x = this->x;
    double y = this->y;
    double z = this->z;
    double ax = axis.x;
    double ay = axis.y;
    double az = axis.z;

    // Compute the rotation matrix
    double m11 = c + ax * ax * t;
    double m12 = ax * ay * t - az * s;
    double m13 = ax * az * t + ay * s;
    double m21 = ay * ax * t + az * s;
    double m22 = c + ay * ay * t;
    double m23 = ay * az * t - ax * s;
    double m31 = az * ax * t - ay * s;
    double m32 = az * ay * t + ax * s;
    double m33 = c + az * az * t;

    // Apply the rotation matrix to the vector
    double n_x = m11 * x + m12 * y + m13 * z;
    double n_y = m21 * x + m22 * y + m23 * z;
    double n_z = m31 * x + m32 * y + m33 * z;

    return Vector3D(n_x, n_y, n_z);
  }
  /// Pretty-print
  std::string __str__() const override
  {
    return "(" + str(x) + ", " + str(y) + ", " + str(z) + ")";
  }
};

struct Vector3DHash
{
  std::size_t operator()(const Vector3D &v) const
  {
    std::hash<double> hasher;
    std::size_t x_hash = hasher(v.x);
    std::size_t y_hash = hasher(v.y);
    std::size_t z_hash = hasher(v.z);
    return x_hash ^ (y_hash << 1) ^ (z_hash << 2);
  }
};

} // namespace DTCC_BUILDER

#endif
