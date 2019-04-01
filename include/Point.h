// Point classes for 2D and 3D.
// Copyright (C) 2018 Anders Logg.

#ifndef VC_POINT_H
#define VC_POINT_H

#include <iostream>
#include <iomanip>
#include <cmath>

#include "Parameters.h"

namespace VirtualCity
{

class Point2D
{
public:

    double x;
    double y;

    Point2D() : x(0), y(0) {}
    Point2D(double x, double y) : x(x), y(y) {}

    Point2D operator+ (const Point2D& p) const
    {
        Point2D q(x + p.x, y + p.y);
        return q;
    }

    Point2D operator+= (const Point2D& p)
    {
        x += p.x;
        y += p.y;
        return *this;
    }

    Point2D operator- (const Point2D& p) const
    {
        Point2D q(x - p.x, y - p.y);
        return q;
    }

    Point2D operator-= (const Point2D& p)
    {
        x -= p.x;
        y -= p.y;
        return *this;
    }

    Point2D operator* (double a) const
    {
        Point2D q(x * a, y * a);
        return q;
    }

    Point2D operator*= (double a)
    {
        x *= a;
        y *= a;
        return *this;
    }

    Point2D operator/ (double a) const
    {
        Point2D q(x / a, y / a);
        return q;
    }

    Point2D operator/= (double a)
    {
        x /= a;
        y /= a;
        return *this;
    }

    double Magnitude() const
    {
        return sqrt(SquaredMagnitude());
    }

    double SquaredMagnitude() const
    {
        return x * x + y * y;
    }

};

class Point3D
{
public:

    double x;
    double y;
    double z;

    Point3D() : x(0), y(0), z(0) {}
    Point3D(double x, double y, double z) : x(x), y(y), z(z) {}

    Point3D operator+ (const Point3D& p) const
    {
        Point3D q(x + p.x, y + p.y, z + p.z);
        return q;
    }

    Point3D operator+= (const Point3D& p)
    {
        x += p.x;
        y += p.y;
        z += p.z;
        return *this;
    }

    Point3D operator- (const Point3D& p) const
    {
        Point3D q(x - p.x, y - p.y, z - p.z);
        return q;
    }

    Point3D operator-= (const Point3D& p)
    {
        x -= p.x;
        y -= p.y;
        z -= p.z;
        return *this;
    }

    Point3D operator* (double a) const
    {
        Point3D q(x * a, y * a, z * a);
        return q;
    }

    Point3D operator*= (double a)
    {
        x *= a;
        y *= a;
        z *= a;
        return *this;
    }

    Point3D operator/ (double a) const
    {
        Point3D q(x / a, y / a, z / a);
        return q;
    }

    Point3D operator/= (double a)
    {
        x /= a;
        y /= a;
        z /= a;
        return *this;
    }

    double Magnitude() const
    {
        return sqrt(SquaredMagnitude());
    }

    double SquaredMagnitude() const
    {
        return x * x + y * y + z * z;
    }

};

std::ostream& operator<<(std::ostream& stream, const Point2D& p)
{
    stream << std::setprecision(Parameters::Precision)
           << "(" << p.x << ", " << p.y << ")";
}

std::ostream& operator<<(std::ostream& stream, const Point3D& p)
{
    stream << std::setprecision(Parameters::Precision)
           << "(" << p.x << ", " << p.y << ", " << p.z << ")";
}

}

#endif
