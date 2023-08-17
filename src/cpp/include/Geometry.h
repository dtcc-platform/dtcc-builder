// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_GEOMETRY_H
#define DTCC_GEOMETRY_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <stack>
#include <tuple>
#include <vector>
#include <iso646.h> 

#include "BoundingBox.h"
#include "Constants.h"
#include "model/Mesh.h"
#include "model/Point.h"
#include "model/Polygon.h"
#include "model/Simplices.h"
#include "model/Vector.h"
#include "model/VolumeMesh.h"

namespace DTCC_BUILDER
{

class Geometry
{
public:
  // FIXME: This needs to be reworked in light of the introduction of two
  // different classes Point and Vector. Currently somewhat inconsistent.

  // Compute squared norm (2D)
  static double SquaredNorm2D(const Vector2D &v) { return Dot2D(v, v); }

  // Compute squared norm (3D)
  static double SquaredNorm3D(const Vector3D &v) { return Dot3D(v, v); }

  // Compute norm (2D)
  static double Norm2D(const Vector2D &v)
  {
    return std::sqrt(SquaredNorm2D(v));
  }

  // Compute norm (3D)
  static double Norm3D(const Vector3D &v)
  {
    return std::sqrt(SquaredNorm3D(v));
  }

  // Compute dot product (2D)
  static double Dot2D(const Vector2D &u, const Vector2D &v)
  {
    return u.x * v.x + u.y * v.y;
  }

  // Compute dot product (3D)
  static double Dot3D(const Vector3D &u, const Vector3D &v)
  {
    return u.x * v.x + u.y * v.y + u.z * v.z;
  }

  // Compute cross product (3D)
  static Vector3D Cross3D(const Vector3D &u, const Vector3D &v)
  {
    return Vector3D(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z,
                    u.x * v.y - u.y * v.x);
  }

  // Compute distance between points (2D)
  static double Distance2D(const Point2D &p, const Point2D &q)
  {
    return std::sqrt(SquaredDistance2D(p, q));
  }

  // Compute distance between segment (p0, p1) and point q (2D)
  static double
  Distance2D(const Point2D &p0, const Point2D &p1, const Point2D &q)
  {
    return std::sqrt(SquaredDistance2D(p0, p1, q));
  }

  // Compute distance between polygon and point (2D)
  static double Distance2D(const Polygon &polygon, const Point2D &p)
  {
    return std::sqrt(SquaredDistance2D(polygon, p));
  }

  // Compute distance between polygons (2D)
  static double Distance2D(const Polygon &polygon0, const Polygon &polygon1)
  {
    return std::sqrt(SquaredDistance2D(polygon0, polygon1));
  }

  // Compute distance between points (3D)
  static double Distance3D(const Point3D &p, const Point3D &q)
  {
    return std::sqrt(SquaredDistance3D(p, q));
  }

  // Compute squared distance between points (2D)
  static double SquaredDistance2D(const Point2D &p, const Point2D &q)
  {
    const double dx = p.x - q.x;
    const double dy = p.y - q.y;
    return dx * dx + dy * dy;
  }

  // Compute squared distance between segment (p0, p1) and point q (2D)
  static double
  SquaredDistance2D(const Point2D &p0, const Point2D &p1, const Point2D &q)
  {
    // Project point to line
    const Vector2D u(p0, q);
    const Vector2D v(p0, p1);
    const Point2D p = p0 + v * (Dot2D(u, v) / v.SquaredMagnitude());

    // Check whether projected point is inside segment. Check either
    // x or y coordinates depending on which is largest (most stable)
    const bool inside =
        std::abs(v.x) > std::abs(v.y)
            ? std::min(p0.x, p1.x) <= p.x && p.x <= std::max(p0.x, p1.x)
            : std::min(p0.y, p1.y) <= p.y && p.y <= std::max(p0.y, p1.y);

    // Use distance to projection if inside
    if (inside)
      return SquaredDistance2D(p, q);

    // Otherwise use distance to closest end point
    const double d0 = SquaredDistance2D(p0, q);
    const double d1 = SquaredDistance2D(p1, q);
    return std::min(d0, d1);
  }

  // Compute squared distance between polygon and point(2D)
  static double SquaredDistance2D(const Polygon &polygon, const Point2D &p)
  {
    // Check if point is contained in polygon
    if (PolygonContains2D(polygon, p))
      return 0.0;

    // If not, compute minimal squared distance to all segments
    double d2Min = std::numeric_limits<double>::max();
    for (size_t i = 0; i < polygon.Vertices.size(); i++)
    {
      Point2D p0 = polygon.Vertices[i];
      Point2D p1 = polygon.Vertices[(i + 1) % polygon.Vertices.size()];
      d2Min = std::min(d2Min, SquaredDistance2D(p0, p1, p));
    }

    return d2Min;
  }

  // Compute squared distance between polygons (2D)
  static double SquaredDistance2D(const Polygon &polygon0,
                                  const Polygon &polygon1)
  {
    //    std::cout << "Checking polygon with " << polygon0.Vertices.size()
    //              << " vertices" << std::endl;
    //    std::cout << "Checking polygon with " << polygon1.Vertices.size()
    //             << " vertices" << std::endl;
    double d2Min = std::numeric_limits<double>::max();

    // Check all vertices in first polygon
    for (auto const &p : polygon0.Vertices)
    {
      //     std::cout << "Point of polygon0 " << str(p) << std::endl;
      d2Min = std::min(d2Min, SquaredDistance2D(polygon1, p));
      //     std::cout << "D2min " << d2Min << std::endl;
    }
    // Check all vertices in second polygon
    for (auto const &p : polygon1.Vertices)
    {
      //      std::cout << "Point of polygon1 " << str(p) << std::endl;
      d2Min = std::min(d2Min, SquaredDistance2D(polygon0, p));
      //      std::cout << "D2min " << d2Min << std::endl;
    }
    return d2Min;
  }

  // Compute squared distance between points (3D)
  static double SquaredDistance3D(const Point3D &p, const Point3D &q)
  {
    const double dx = p.x - q.x;
    const double dy = p.y - q.y;
    const double dz = p.z - q.z;
    return dx * dx + dy * dy + dz * dz;
  }

  // Compute orientation of point q relative to edge (p0, p1) (2D)
  static double Orient2D(const Point2D &p0, const Point2D &p1, const Point2D &q)
  {
    const Vector2D u(p0, p1);
    const Vector2D v(p0, q);
    return u.x * v.y - u.y * v.x;
  }

  // Compute sign of point q relative to edge (p0, p1) (2D).
  // -1 --- (p0) --- 0 --- (p1) --- +1
  static int EdgeSign2D(const Point2D &p0, const Point2D &p1, const Point2D &q)
  {
    double l{}, d0{}, d1{};
    if (std::abs(p0.x - p1.x) > std::abs(p0.y - p1.y))
    {
      l = std::abs(p0.x - p1.x);
      d0 = std::abs(p0.x - q.x);
      d1 = std::abs(p1.x - q.x);
    }
    else
    {
      l = std::abs(p0.y - p1.y);
      d0 = std::abs(p0.y - q.y);
      d1 = std::abs(p1.y - q.y);
    }
    if (d0 > l - Constants::Epsilon && d0 > d1)
      return 1;
    else if (d1 > l - Constants::Epsilon && d1 > d0)
      return -1;
    else
      return 0;
  }

  // Compute strictly increasing function [-pi, pi] -> [-2, 2] of angle of v
  // relative to u. This is a cheap alternative compared to working with asin,
  // acos.
  static double VectorAngle2D(const Vector2D &u, const Vector2D &v)
  {
    const double u2 = u.x * u.x + u.y * u.y;
    const double v2 = v.x * v.x + v.y * v.y;
    const double sin = u.x * v.y - u.y * v.x;
    const double cos = u.x * v.x + u.y * v.y;
    const double a = sin * sin / (u2 * v2);
    return (sin > 0.0 ? (cos > 0.0 ? a : 2.0 - a) : (cos > 0.0 ? -a : a - 2.0));
  }

  // Compute quadrant angle of point p relative to polygon (2D)
  static int QuadrantAngle2D(const Point2D &p,
                             const std::vector<Point2D> &polygon)
  {
    // Compute angle to first vertex
    Point2D q0 = polygon[0];
    int v0 = QuadrantAngle2D(q0, p);

    // Sum up total angle
    int totalAngle = 0;
    for (size_t i = 1; i < polygon.size() + 1; i++)
    {
      // Compute angle increment
      Point2D q1 = polygon[i % polygon.size()];
      int v1 = QuadrantAngle2D(q1, p);
      int dv = v1 - v0;

      // Adjust angle increment for wrap-around
      if (dv == 3)
        dv = -1;
      else if (dv == -3)
        dv = 1;
      else if (dv == 2 || dv == -2)
      {
        double xx = q1.x - ((q1.y - p.y) * ((q0.x - q1.x) / (q0.y - q1.y)));
        if (xx > p.x)
          dv = -dv;
      }

      // Add to total angle and update
      totalAngle += dv;
      q0 = q1;
      v0 = v1;
    }

    return totalAngle;
  }

  // Compute quadrant angle of point p relative to point q (2D)
  static int QuadrantAngle2D(const Point2D &p, const Point2D &q)
  {
    return ((p.x > q.x) ? ((p.y > q.y) ? 0 : 3) : ((p.y > q.y) ? 1 : 2));
  }

  // Compute face normal
  static Vector3D FaceNormal3D(const Simplex2D &face, const VolumeMesh &mesh3D)
  {
    const Vector3D p0{mesh3D.Vertices[face.v0]};
    const Vector3D p1{mesh3D.Vertices[face.v1]};
    const Vector3D p2{mesh3D.Vertices[face.v2]};
    const Vector3D u = p1 - p0;
    const Vector3D v = p2 - p0;
    Vector3D n = Cross3D(u, v);
    n /= Geometry::Norm3D(n);
    return n;
  }

  // Compute cell center
  static Point3D CellCenter3D(const Simplex3D &cell, const VolumeMesh &mesh3D)
  {
    Vector3D c{};
    c += Vector3D(mesh3D.Vertices[cell.v0]);
    c += Vector3D(mesh3D.Vertices[cell.v1]);
    c += Vector3D(mesh3D.Vertices[cell.v2]);
    c += Vector3D(mesh3D.Vertices[cell.v3]);
    c /= 4.0;
    return c;
  }

  // Compute signed determinant of polygon (2D)
  static double PolygonDeterminant2D(const Polygon &polygon)
  {
    double sum = 0.0;
    for (size_t i = 0; i < polygon.Vertices.size(); i++)
    {
      Point2D p0 = polygon.Vertices[i];
      Point2D p1 = polygon.Vertices[(i + 1) % polygon.Vertices.size()];
      sum += (p1.x - p0.x) * (p1.y + p0.y);
    }
    return sum;
  }

  // Compute orientation of polygon (0 = counter-clockwise, 1 = clockwise)
  static size_t PolygonOrientation2D(const Polygon &polygon)
  {
    return PolygonDeterminant2D(polygon) < 0 ? 0 : 1;
  }

  // Compute area of polygon (2D)
  static double PolygonArea(const Polygon &polygon)
  {
    return 0.5 * std::abs(PolygonDeterminant2D(polygon));
  }

  // Computer Perimeter of polygon (2D)
  static double PolygonPerimeter2D(const Polygon &polygon)
  {
    double sum = 0.0;
    for (size_t i = 0; i < polygon.Vertices.size(); i++)
    {
      Point2D p0 = polygon.Vertices[i];
      Point2D p1 = polygon.Vertices[(i + 1) % polygon.Vertices.size()];
      sum += Geometry::Distance2D(p0, p1);
    }
    return sum;
  }

  // Compute center of polygon (2D)
  static Point2D PolygonCenter2D(const Polygon &polygon)
  {
    Vector2D o{};
    Vector2D c{};
    for (auto const &p : polygon.Vertices)
      c += Vector2D(o, p);
    c /= static_cast<double>(polygon.Vertices.size());
    return c;
  }

  // Compute radius of polygon relative to center (2D)
  static double PolygonRadius2D(const Polygon &polygon, const Point2D &center)
  {
    double r2max = 0.0;
    for (auto const &p : polygon.Vertices)
    {
      const double r2 = SquaredDistance2D(p, center);
      if (r2 > r2max)
        r2max = r2;
    }
    return std::sqrt(r2max);
  }

  // FIXME: Cleanup Contains vs Intersects vs Collide. What should we name it?

  // Check whether edge (p0, p1) contains point q. It is assumed that the
  // point is located on the line defined by the edge.
  static bool EdgeContains2D(const Point2D &p0,
                             const Point2D &p1,
                             const Point2D &q,
                             double tol = 0.0)
  {
    const Vector2D v(p0, p1);
    if (std::abs(v.x) > std::abs(v.y))
      return std::min(p0.x, p1.x) - tol < q.x and
             std::max(p0.x, p1.x) + tol > q.x;
    else
      return std::min(p0.y, p1.y) - tol < q.y and
             std::max(p0.y, p1.y) + tol > q.y;
  }

  // Check whether polygon contains point (2D)
  static bool PolygonContains2D(const Polygon &polygon, const Point2D &p)
  {
    // Compute total quadrant relative to polygon. If the point
    // is inside the polygon, the angle should be 4 (or -4).
    return Geometry::QuadrantAngle2D(p, polygon.Vertices) != 0;
  }

  // Check whether bounding box contains point (2D)
  static bool BoundingBoxContains2D(const BoundingBox2D &bbox,
                                    const Point2D &p,
                                    double margin = 0.0)
  {
    return (bbox.P.x + margin <= p.x && p.x + margin <= bbox.Q.x &&
            bbox.P.y + margin <= p.y && p.y + margin <= bbox.Q.y);
  }

  // Check whether bounding box contains point (3D)
  static bool BoundingBoxContains3D(const BoundingBox3D &bbox,
                                    const Point3D &p,
                                    double margin = 0.0)
  {
    return (bbox.P.x + margin <= p.x && p.x + margin <= bbox.Q.x &&
            bbox.P.y + margin <= p.y && p.y + margin <= bbox.Q.y &&
            bbox.P.z + margin <= p.z && p.z + margin <= bbox.Q.z);
  }

  // Check whether bounding box contains polygon (2D)
  static bool BoundingBoxContains2D(const BoundingBox2D &bbox,
                                    const Polygon &polygon,
                                    double margin = 0.0)
  {
    for (const auto &p : polygon.Vertices)
      if (!BoundingBoxContains2D(bbox, p, margin))
        return false;
    return true;
  }

  // Check whether edges (p0, p1) and (q0, q1) intersect
  static bool Intersects2D(const Point2D &p0,
                           const Point2D &p1,
                           const Point2D &q0,
                           const Point2D &q1)
  {
    return (Orient2D(p0, p1, q0) * Orient2D(p0, p1, q1) <= 0.0 &&
            Orient2D(q0, q1, p0) * Orient2D(q0, q1, p1) <= 0.0);
  }

  // Check wheter bounding boxes intersect (2D)
  static bool Intersect2D(const BoundingBox2D &bboxA,
                          const BoundingBox2D &bboxB)
  {
    return (bboxA.P.x <= bboxB.Q.x && bboxB.P.x <= bboxA.Q.x &&
            bboxA.P.y <= bboxB.Q.y && bboxB.P.y <= bboxA.Q.y);
  }

  // Check whether polygon intersects with polygon (2D)
  static bool Intersects2D(const Polygon &polygonA, const Polygon &polygonB)
  {
    // Check if bounding boxes intersect
    if (!Intersect2D(BoundingBox2D(polygonA), BoundingBox2D(polygonB)))
      return false;

    // Check if any edge of polygonA intersects with any edge of polygonB
    for (const auto &p0 : polygonA.Vertices)
      for (const auto &p1 : polygonA.Vertices)
        for (const auto &q0 : polygonB.Vertices)
          for (const auto &q1 : polygonB.Vertices)
            if (Intersects2D(p0, p1, q0, q1))
              return true;

    // Check if one polygon contains the other.
    if (PolygonContains2D(polygonA, polygonB.Vertices[0]) ||
        PolygonContains2D(polygonB, polygonA.Vertices[0]))
      return true;

    return false;
  }

  // Compute intersection between edges p0 - p1 and q0 - q1 (2D)
  static Point2D EdgeIntersection2D(const Point2D &p0,
                                    const Point2D &p1,
                                    const Point2D &q0,
                                    const Point2D &q1)
  {
    // Solve for intersection: p0 + k*(p1 - p0) = q0 + l*(q1 - q0)

    // Compute vectors
    const Vector2D u(p0, q0);
    const Vector2D v(p0, p1);
    const Vector2D w(q0, q1);

    // Create linear system
    const double a = v.x;
    const double b = -w.x;
    const double c = v.y;
    const double d = -w.y;
    const double e = u.x;
    const double f = u.y;

    // Compute determinant
    const double det = a * d - b * c;

    // Check if close to parallel
    // if (std::(det) < parallelTolerance)
    // {
    //     throw std::runtime_error("Segments are parallel.");
    // }

    // Solve linear system
    const double k = (d * e - b * f) / det;
    Point2D p = p0 + v * k;

    return p;
  }

  // Compute convex hull of point set (2D)
  static Polygon ConvexHull2D(const std::vector<Point2D> &points)
  {
    // The convex hull is computed by doing a Graham scan: select an
    // extreme base point, sort remaining points by angle and then
    // add points that create a left turn around the perimeter.

    // Find point with smallest y-coordinate. If y-coordinate is
    // the same, sort by smallest x-coordinate.
    double xMin = points[0].x;
    double yMin = points[0].y;
    size_t iMin = 0;
    const size_t numPoints = points.size();
    for (size_t i = 1; i < numPoints; i++)
    {
      const double x = points[i].x;
      const double y = points[i].y;
      if (y < yMin || (y == yMin && x < xMin))
      {
        xMin = x;
        yMin = y;
        iMin = i;
      }
    }

    // Set base point
    const size_t baseIndex = iMin;
    Point2D basePoint = points[baseIndex];

    // Compute angles and distances relative to base point
    std::vector<std::tuple<double, double, size_t>> angles(numPoints - 1);
    size_t k = 0;
    for (size_t i = 0; i < numPoints; i++)
    {
      // Skip base point
      if (i == baseIndex)
        continue;

      // Compute angle (negative cosine) and distance
      const Point2D &p = points[i];
      const Vector2D v(basePoint, p);
      const double distance = v.Magnitude();
      const double angle =
          (distance > Constants::Epsilon ? -v.x / distance : 0.0);

      // Store angle and distance along with index (for sorting)
      angles[k++] = std::make_tuple(angle, distance, i);
    }

    // Sort by angles (primary) and distance (secondary) to base point
    std::sort(angles.begin(), angles.end());

    // Filter out points with unique angles, keeping only furthest point
    std::vector<size_t> filteredIndices;
    double lastAngle = 2.0; // no angle has this value
    for (size_t i = 0; i < numPoints - 1; i++)
    {
      // Get data for current point
      const double currentAngle = std::get<0>(angles[i]);
      const size_t currentIndex = std::get<2>(angles[i]);

      // Add point or replace last point
      if (std::abs(currentAngle - lastAngle) > Constants::Epsilon)
        filteredIndices.push_back(currentIndex);
      else
        filteredIndices[filteredIndices.size() - 1] = currentIndex;

      // Update last index
      lastAngle = currentAngle;
    }

    // Create stack of points and add first three candidates
    std::stack<size_t> convexHull;
    convexHull.push(baseIndex);
    convexHull.push(filteredIndices[0]);
    convexHull.push(filteredIndices[1]);

    // Graham-Scan: Push candidates to stack and pop until
    // we have a left turn
    for (size_t i = 2; i < filteredIndices.size(); i++)
    {
      // Get next point
      const size_t i2 = filteredIndices[i];
      const Point2D &p2 = points[i2];

      // Keep popping from stack until we see a left turn
      while (true)
      {
        // Get last two points from stack
        const size_t i1 = convexHull.top();
        convexHull.pop();
        const size_t i0 = convexHull.top();
        const Point2D &p0 = points[i0];
        const Point2D &p1 = points[i1];

        // Check orientation, keep p1 if orientation is positive
        if (Orient2D(p0, p1, p2) > Constants::Epsilon)
        {
          convexHull.push(i1);
          break;
        }
      }

      // Push next candidate to stack
      convexHull.push(i2);
    }

    // Extract polygon points from stack
    Polygon polygon;
    while (!convexHull.empty())
    {
      polygon.Vertices.push_back(points[convexHull.top()]);
      convexHull.pop();
    }

    // Reverse polygon to make it counter-clockwise
    std::reverse(polygon.Vertices.begin(), polygon.Vertices.end());

    return polygon;
  }



};

} // namespace DTCC_BUILDER

#endif
