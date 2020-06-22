// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_VECTOR_H
#define DTCC_VECTOR_H

#include <cmath>

#include "Logging.h"

// FIXME: Temporary conversion during transition to Point/Vector
#include "Point.h"

namespace DTCC
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
    Vector2D() {}

    /// Create vector with given components.
    ///
    /// @param x First component
    /// @param y Second component
    Vector2D(double x, double y) : x(x), y(y) {}

    /// Create vector between points.
    ///
    /// @param p First point
    /// @param q Second point
    Vector2D(const Point2D &p, const Point2D &q) : x(q.x - p.x), y(q.y - p.y) {}

    // FIXME: Temporary conversion during transition to Point/Vector
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
      Vector2D q(x - p.x, y - p.y);
      return q;
    }

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

    double Magnitude() const { return sqrt(SquaredMagnitude()); }

    double SquaredMagnitude() const { return x * x + y * y; }

    void Normalize() { (*this) /= Magnitude(); }

    /// Pretty-print
    std::string __str__() const
    {
      return "<" + str(x) + ", " + str(y) + ">";
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
    Vector3D() {}

    /// Create vector with given components.
    ///
    /// @param x First component
    /// @param y Second component
    /// @param z Third component
    Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}

    /// Create vector between points.
    ///
    /// @param p First point
    /// @param q Second point
    Vector3D(const Point3D &p, const Point3D &q)
        : x(q.x - p.x), y(q.y - p.y), z(q.z - p.z)
    {
    }

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
      Vector3D q(x - p.x, y - p.y, z - p.z);
      return q;
    }

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

    double Magnitude() const { return sqrt(SquaredMagnitude()); }

    double SquaredMagnitude() const { return x * x + y * y + z * z; }

    /// Pretty-print
    std::string __str__() const
    {
      return "<" + str(x) + ", " + str(y) + ", " + str(z) + ">";
    }

  };

} // namespace DTCC

#endif
