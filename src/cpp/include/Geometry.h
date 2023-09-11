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
  static double squared_norm_2d(const Vector2D &v) { return dot_2d(v, v); }

  // Compute squared norm (3D)
  static double squared_norm_3d(const Vector3D &v) { return dot_3d(v, v); }

  // Compute norm (2D)
  static double norm_2d(const Vector2D &v)
  {
    return std::sqrt(squared_norm_2d(v));
  }

  // Compute norm (3D)
  static double norm_3d(const Vector3D &v)
  {
    return std::sqrt(squared_norm_3d(v));
  }

  // Compute dot product (2D)
  static double dot_2d(const Vector2D &u, const Vector2D &v)
  {
    return u.x * v.x + u.y * v.y;
  }

  // Compute dot product (3D)
  static double dot_3d(const Vector3D &u, const Vector3D &v)
  {
    return u.x * v.x + u.y * v.y + u.z * v.z;
  }

  // Compute cross product (3D)
  static Vector3D cross_3d(const Vector3D &u, const Vector3D &v)
  {
    return Vector3D(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z,
                    u.x * v.y - u.y * v.x);
  }

  // Compute distance between points (2D)
  static double distance_2d(const Vector2D &p, const Vector2D &q)
  {
    return std::sqrt(squared_distance_2d(p, q));
  }

  // Compute distance between segment (p0, p1) and point q (2D)
  static double
  distance_2d(const Vector2D &p0, const Vector2D &p1, const Vector2D &q)
  {
    return std::sqrt(squared_distance_2d(p0, p1, q));
  }

  // Compute distance between polygon and point (2D)
  static double distance_2d(const Polygon &polygon, const Vector2D &p)
  {
    return std::sqrt(squared_distance_2d(polygon, p));
  }

  // Compute distance between polygons (2D)
  static double distance_2d(const Polygon &polygon0, const Polygon &polygon1)
  {
    return std::sqrt(squared_distance_2d(polygon0, polygon1));
  }

  // Compute distance between points (3D)
  static double distance_3d(const Vector3D &p, const Vector3D &q)
  {
    return std::sqrt(squared_distance_3d(p, q));
  }

  // Compute squared distance between points (2D)
  static double squared_distance_2d(const Vector2D &p, const Vector2D &q)
  {
    const double dx = p.x - q.x;
    const double dy = p.y - q.y;
    return dx * dx + dy * dy;
  }

  // Compute squared distance between segment (p0, p1) and point q (2D)
  static double
  squared_distance_2d(const Vector2D &p0, const Vector2D &p1, const Vector2D &q)
  {
    // Project point to line
    const Vector2D u(p0, q);
    const Vector2D v(p0, p1);
    const Vector2D p = p0 + v * (dot_2d(u, v) / v.squared_magnitude());

    // Check whether projected point is inside segment. Check either
    // x or y coordinates depending on which is largest (most stable)
    const bool inside =
        std::abs(v.x) > std::abs(v.y)
            ? std::min(p0.x, p1.x) <= p.x && p.x <= std::max(p0.x, p1.x)
            : std::min(p0.y, p1.y) <= p.y && p.y <= std::max(p0.y, p1.y);

    // Use distance to projection if inside
    if (inside)
      return squared_distance_2d(p, q);

    // Otherwise use distance to closest end point
    const double d0 = squared_distance_2d(p0, q);
    const double d1 = squared_distance_2d(p1, q);
    return std::min(d0, d1);
  }

  // Compute squared distance between polygon and point(2D)
  static double squared_distance_2d(const Polygon &polygon, const Vector2D &p)
  {
    // Check if point is contained in polygon
    if (polygon_contains_2d(polygon, p))
      return 0.0;

    // If not, compute minimal squared distance to all segments
    double d_to_min = std::numeric_limits<double>::max();
    for (size_t i = 0; i < polygon.vertices.size(); i++)
    {
      Vector2D p0 = polygon.vertices[i];
      Vector2D p1 = polygon.vertices[(i + 1) % polygon.vertices.size()];
      d_to_min = std::min(d_to_min, squared_distance_2d(p0, p1, p));
    }

    return d_to_min;
  }

  // Compute squared distance between polygons (2D)
  static double squared_distance_2d(const Polygon &polygon0,
                                    const Polygon &polygon1)
  {
    //    std::cout << "Checking polygon with " << polygon0.vertices.size()
    //              << " vertices" << std::endl;
    //    std::cout << "Checking polygon with " << polygon1.vertices.size()
    //             << " vertices" << std::endl;
    double d_to_min = std::numeric_limits<double>::max();

    // Check all vertices in first polygon
    for (auto const &p : polygon0.vertices)
    {
      //     std::cout << "Point of polygon0 " << str(p) << std::endl;
      d_to_min = std::min(d_to_min, squared_distance_2d(polygon1, p));
      //     std::cout << "D2min " << d_to_min << std::endl;
    }
    // Check all vertices in second polygon
    for (auto const &p : polygon1.vertices)
    {
      //      std::cout << "Point of polygon1 " << str(p) << std::endl;
      d_to_min = std::min(d_to_min, squared_distance_2d(polygon0, p));
      //      std::cout << "D2min " << d_to_min << std::endl;
    }
    return d_to_min;
  }

  // Compute squared distance between points (3D)
  static double squared_distance_3d(const Vector3D &p, const Vector3D &q)
  {
    const double dx = p.x - q.x;
    const double dy = p.y - q.y;
    const double dz = p.z - q.z;
    return dx * dx + dy * dy + dz * dz;
  }

  // Compute orientation of point q relative to edge (p0, p1) (2D)
  static double
  orient_2d(const Vector2D &p0, const Vector2D &p1, const Vector2D &q)
  {
    const Vector2D u(p0, p1);
    const Vector2D v(p0, q);
    return u.x * v.y - u.y * v.x;
  }

  // Compute sign of point q relative to edge (p0, p1) (2D).
  // -1 --- (p0) --- 0 --- (p1) --- +1
  static int
  edge_sign_2d(const Vector2D &p0, const Vector2D &p1, const Vector2D &q)
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
    if (d0 > l - Constants::epsilon && d0 > d1)
      return 1;
    else if (d1 > l - Constants::epsilon && d1 > d0)
      return -1;
    else
      return 0;
  }

  // Compute strictly increasing function [-pi, pi] -> [-2, 2] of angle of v
  // relative to u. This is a cheap alternative compared to working with asin,
  // acos.
  static double vector_angle_2d(const Vector2D &u, const Vector2D &v)
  {
    const double u2 = u.x * u.x + u.y * u.y;
    const double v2 = v.x * v.x + v.y * v.y;
    const double sin = u.x * v.y - u.y * v.x;
    const double cos = u.x * v.x + u.y * v.y;
    const double a = sin * sin / (u2 * v2);
    return (sin > 0.0 ? (cos > 0.0 ? a : 2.0 - a) : (cos > 0.0 ? -a : a - 2.0));
  }

  // Compute quadrant angle of point p relative to polygon (2D)
  static int quadrant_angle_2d(const Vector2D &p,
                               const std::vector<Vector2D> &polygon)
  {
    // Compute angle to first vertex
    Vector2D q0 = polygon[0];
    int v0 = quadrant_angle_2d(q0, p);

    // Sum up total angle
    int total_angle = 0;
    for (size_t i = 1; i < polygon.size() + 1; i++)
    {
      // Compute angle increment
      Vector2D q1 = polygon[i % polygon.size()];
      int v1 = quadrant_angle_2d(q1, p);
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
      total_angle += dv;
      q0 = q1;
      v0 = v1;
    }

    return total_angle;
  }

  // Compute quadrant angle of point p relative to point q (2D)
  static int quadrant_angle_2d(const Vector2D &p, const Vector2D &q)
  {
    return ((p.x > q.x) ? ((p.y > q.y) ? 0 : 3) : ((p.y > q.y) ? 1 : 2));
  }

  // Compute face normal
  static Vector3D face_normal_3d(const Simplex2D &face,
                                 const VolumeMesh &mesh_3d)
  {
    const Vector3D p0{mesh_3d.vertices[face.v0]};
    const Vector3D p1{mesh_3d.vertices[face.v1]};
    const Vector3D p2{mesh_3d.vertices[face.v2]};
    const Vector3D u = p1 - p0;
    const Vector3D v = p2 - p0;
    Vector3D n = cross_3d(u, v);
    n /= Geometry::norm_3d(n);
    return n;
  }

  // Compute cell center
  static Vector3D cell_center_3d(const Simplex3D &cell,
                                 const VolumeMesh &mesh_3d)
  {
    Vector3D c{};
    c += Vector3D(mesh_3d.vertices[cell.v0]);
    c += Vector3D(mesh_3d.vertices[cell.v1]);
    c += Vector3D(mesh_3d.vertices[cell.v2]);
    c += Vector3D(mesh_3d.vertices[cell.v3]);
    c /= 4.0;
    return c;
  }

  // Compute signed determinant of polygon (2D)
  static double polygon_determinant_2d(const Polygon &polygon)
  {
    double sum = 0.0;
    for (size_t i = 0; i < polygon.vertices.size(); i++)
    {
      Vector2D p0 = polygon.vertices[i];
      Vector2D p1 = polygon.vertices[(i + 1) % polygon.vertices.size()];
      sum += (p1.x - p0.x) * (p1.y + p0.y);
    }
    return sum;
  }

  // Compute orientation of polygon (0 = counter-clockwise, 1 = clockwise)
  static size_t polygon_orientation_2d(const Polygon &polygon)
  {
    return polygon_determinant_2d(polygon) < 0 ? 0 : 1;
  }

  // Compute area of polygon (2D)
  static double polygon_area(const Polygon &polygon)
  {
    return 0.5 * std::abs(polygon_determinant_2d(polygon));
  }

  // Computer Perimeter of polygon (2D)
  static double polygon_perimeter_2d(const Polygon &polygon)
  {
    double sum = 0.0;
    for (size_t i = 0; i < polygon.vertices.size(); i++)
    {
      Vector2D p0 = polygon.vertices[i];
      Vector2D p1 = polygon.vertices[(i + 1) % polygon.vertices.size()];
      sum += Geometry::distance_2d(p0, p1);
    }
    return sum;
  }

  // Compute center of polygon (2D)
  static Vector2D polygon_center_2d(const Polygon &polygon)
  {
    Vector2D o{};
    Vector2D c{};
    for (auto const &p : polygon.vertices)
      c += Vector2D(o, p);
    c /= static_cast<double>(polygon.vertices.size());
    return c;
  }

  // Compute radius of polygon relative to center (2D)
  static double polygon_radius_2d(const Polygon &polygon,
                                  const Vector2D &center)
  {
    double r_to_max = 0.0;
    for (auto const &p : polygon.vertices)
    {
      const double r2 = squared_distance_2d(p, center);
      if (r2 > r_to_max)
        r_to_max = r2;
    }
    return std::sqrt(r_to_max);
  }

  // FIXME: Cleanup Contains vs Intersects vs Collide. What should we name it?

  // Check whether edge (p0, p1) contains point q. It is assumed that the
  // point is located on the line defined by the edge.
  static bool edge_contains_2d(const Vector2D &p0,
                               const Vector2D &p1,
                               const Vector2D &q,
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
  static bool polygon_contains_2d(const Polygon &polygon, const Vector2D &p)
  {
    // Compute total quadrant relative to polygon. If the point
    // is inside the polygon, the angle should be 4 (or -4).
    return Geometry::quadrant_angle_2d(p, polygon.vertices) != 0;
  }

  // Check whether bounding box contains point (2D)
  static bool bounding_box_contains_2d(const BoundingBox2D &bbox,
                                       const Vector2D &p,
                                       double margin = 0.0)
  {
    return (bbox.P.x + margin <= p.x && p.x + margin <= bbox.Q.x &&
            bbox.P.y + margin <= p.y && p.y + margin <= bbox.Q.y);
  }

  // Check whether bounding box contains point (3D)
  static bool bounding_box_contains_3d(const BoundingBox3D &bbox,
                                       const Vector3D &p,
                                       double margin = 0.0)
  {
    return (bbox.P.x + margin <= p.x && p.x + margin <= bbox.Q.x &&
            bbox.P.y + margin <= p.y && p.y + margin <= bbox.Q.y &&
            bbox.P.z + margin <= p.z && p.z + margin <= bbox.Q.z);
  }

  // Check whether bounding box contains polygon (2D)
  static bool bounding_box_contains_2d(const BoundingBox2D &bbox,
                                       const Polygon &polygon,
                                       double margin = 0.0)
  {
    for (const auto &p : polygon.vertices)
      if (!bounding_box_contains_2d(bbox, p, margin))
        return false;
    return true;
  }

  // Check whether edges (p0, p1) and (q0, q1) intersect
  static bool intersects_2d(const Vector2D &p0,
                            const Vector2D &p1,
                            const Vector2D &q0,
                            const Vector2D &q1)
  {
    return (orient_2d(p0, p1, q0) * orient_2d(p0, p1, q1) <= 0.0 &&
            orient_2d(q0, q1, p0) * orient_2d(q0, q1, p1) <= 0.0);
  }

  // Check wheter bounding boxes intersect (2D)
  static bool intersect_2d(const BoundingBox2D &bbox_a,
                           const BoundingBox2D &bbox_b)
  {
    return (bbox_a.P.x <= bbox_b.Q.x && bbox_b.P.x <= bbox_a.Q.x &&
            bbox_a.P.y <= bbox_b.Q.y && bbox_b.P.y <= bbox_a.Q.y);
  }

  // Check whether polygon intersects with polygon (2D)
  static bool intersects_2d(const Polygon &polygon_a, const Polygon &polygon_b)
  {
    // Check if bounding boxes intersect
    if (!intersect_2d(BoundingBox2D(polygon_a), BoundingBox2D(polygon_b)))
      return false;

    // Check if any edge of polygon_a intersects with any edge of polygon_b
    for (const auto &p0 : polygon_a.vertices)
      for (const auto &p1 : polygon_a.vertices)
        for (const auto &q0 : polygon_b.vertices)
          for (const auto &q1 : polygon_b.vertices)
            if (intersects_2d(p0, p1, q0, q1))
              return true;

    // Check if one polygon contains the other.
    if (polygon_contains_2d(polygon_a, polygon_b.vertices[0]) ||
        polygon_contains_2d(polygon_b, polygon_a.vertices[0]))
      return true;

    return false;
  }

  // Compute intersection between edges p0 - p1 and q0 - q1 (2D)
  static Vector2D edge_intersection_2d(const Vector2D &p0,
                                       const Vector2D &p1,
                                       const Vector2D &q0,
                                       const Vector2D &q1)
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
    Vector2D p = p0 + v * k;

    return p;
  }

  // Compute convex hull of point set (2D)
  static Polygon convex_hull_2d(const std::vector<Vector2D> &points)
  {
    // The convex hull is computed by doing a Graham scan: select an
    // extreme base point, sort remaining points by angle and then
    // add points that create a left turn around the perimeter.

    // find point with smallest y-coordinate. If y-coordinate is
    // the same, sort by smallest x-coordinate.
    double x_min = points[0].x;
    double y_min = points[0].y;
    size_t i_min = 0;
    const size_t num_points = points.size();
    for (size_t i = 1; i < num_points; i++)
    {
      const double x = points[i].x;
      const double y = points[i].y;
      if (y < y_min || (y == y_min && x < x_min))
      {
        x_min = x;
        y_min = y;
        i_min = i;
      }
    }

    // Set base point
    const size_t base_index = i_min;
    Vector2D base_point = points[base_index];

    // Compute angles and distances relative to base point
    std::vector<std::tuple<double, double, size_t>> angles(num_points - 1);
    size_t k = 0;
    for (size_t i = 0; i < num_points; i++)
    {
      // Skip base point
      if (i == base_index)
        continue;

      // Compute angle (negative cosine) and distance
      const Vector2D &p = points[i];
      const Vector2D v(base_point, p);
      const double distance = v.magnitude();
      const double angle =
          (distance > Constants::epsilon ? -v.x / distance : 0.0);

      // Store angle and distance along with index (for sorting)
      angles[k++] = std::make_tuple(angle, distance, i);
    }

    // Sort by angles (primary) and distance (secondary) to base point
    std::sort(angles.begin(), angles.end());

    // Filter out points with unique angles, keeping only furthest point
    std::vector<size_t> filtered_indices;
    double last_angle = 2.0; // no angle has this value
    for (size_t i = 0; i < num_points - 1; i++)
    {
      // Get data for current point
      const double current_angle = std::get<0>(angles[i]);
      const size_t current_index = std::get<2>(angles[i]);

      // Add point or replace last point
      if (std::abs(current_angle - last_angle) > Constants::epsilon)
        filtered_indices.push_back(current_index);
      else
        filtered_indices[filtered_indices.size() - 1] = current_index;

      // Update last index
      last_angle = current_angle;
    }

    // Create stack of points and add first three candidates
    std::stack<size_t> convex_hull;
    convex_hull.push(base_index);
    convex_hull.push(filtered_indices[0]);
    convex_hull.push(filtered_indices[1]);

    // Graham-Scan: Push candidates to stack and pop until
    // we have a left turn
    for (size_t i = 2; i < filtered_indices.size(); i++)
    {
      // Get next point
      const size_t i2 = filtered_indices[i];
      const Vector2D &p2 = points[i2];

      // Keep popping from stack until we see a left turn
      while (true)
      {
        // Get last two points from stack
        const size_t i1 = convex_hull.top();
        convex_hull.pop();
        const size_t i0 = convex_hull.top();
        const Vector2D &p0 = points[i0];
        const Vector2D &p1 = points[i1];

        // Check orientation, keep p1 if orientation is positive
        if (orient_2d(p0, p1, p2) > Constants::epsilon)
        {
          convex_hull.push(i1);
          break;
        }
      }

      // Push next candidate to stack
      convex_hull.push(i2);
    }

    // Extract polygon points from stack
    Polygon polygon;
    while (!convex_hull.empty())
    {
      polygon.vertices.push_back(points[convex_hull.top()]);
      convex_hull.pop();
    }

    // Reverse polygon to make it counter-clockwise
    std::reverse(polygon.vertices.begin(), polygon.vertices.end());

    return polygon;
  }



};

} // namespace DTCC_BUILDER

#endif
