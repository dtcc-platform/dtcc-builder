// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_POLYFIX_H
#define DTCC_POLYFIX_H

#include "Geometry.h"
#include "Polygon.h"

namespace DTCC
{
/// Polyfix provides algorithms for processing polygons, including
/// polygon cleaning and polygon merging.
class Polyfix
{
public:
  /// Make polygon closed (close polygon when encountering first duplicate
  /// vertex).
  ///
  /// @param polygon The polygon
  /// @param tol Tolerance for small distance
  /// @return 0 if already closed, 1 if modified
  static size_t MakeClosed(Polygon &polygon, double tol)
  {
    // Avoid using sqrt for efficiency
    const double tol2 = tol * tol;

    // Check each vertex against first vertex
    const Vector2D &p0 = polygon.Vertices[0];
    size_t end = 0;
    for (size_t i = 1; i < polygon.Vertices.size() && end == 0; i++)
    {
      // Compute distance to first vertex
      const Vector2D &p = polygon.Vertices[i];
      const double d2 = Geometry::SquaredDistance2D(p, p0);

      // Remove if distance is small
      if (d2 < tol2)
      {
        end = i;
        break;
      }
    }

    // Return if no vertices should be removed
    if (end == 0)
      return 0;

    // Remove vertices
    RemoveVertices(polygon, end);

    return 1;
  }

  /// Make polygon counter-clockwise oriented.
  ///
  /// @param polygon The polygon
  /// @param tol Tolerance
  /// @return 0 if already counter-clockwise, 1 if modified
  static size_t MakeOriented(Polygon &polygon)
  {
    // Return if already counter-clockwise
    if (Geometry::PolygonOrientation2D(polygon) == 0)
      return 0;

    // Reverse polygon
    std::reverse(polygon.Vertices.begin(), polygon.Vertices.end());

    return 1;
  }

  /// Make polygon simple (remove consecutive parallel edges).
  ///
  /// @param polygon The polygon
  /// @param tol Tolerance for small angle (sin of angle)
  /// @return 0 if already simple, 1 if modified
  static bool MakeSimple(Polygon &polygon, double tol)
  {
    // Avoid using sqrt for efficiency
    const double tol2 = tol * tol;

    // Vertices to be removed
    std::vector<size_t> remove;

    // Check each edge
    const size_t numVertices = polygon.Vertices.size();
    for (size_t i = 0; i < numVertices; i++)
    {
      // Get previous, current and next points
      const Vector2D &p0 =
          polygon.Vertices[(i + numVertices - 1) % numVertices];
      const Vector2D &p1 = polygon.Vertices[i];
      const Vector2D &p2 = polygon.Vertices[(i + 1) % numVertices];

      // Compute edges and dot products
      const Vector2D u = p1 - p0;
      const Vector2D v = p2 - p1;
      const double u2 = Geometry::Dot2D(u, u);
      const double v2 = Geometry::Dot2D(v, v);
      const double uv = Geometry::Dot2D(u, v);

      // Remove if angle is small
      if (uv * uv > (1.0 - tol2) * u2 * v2)
        remove.push_back(i);
    }

    // Return if no vertices should be removed
    if (remove.size() == 0)
      return 0;

    // Remove vertices
    RemoveVertices(polygon, remove);

    return 1;
  }

  /// Transform polygon by subtracting given origin.
  ///
  /// @param polygon The polygon
  /// @param origin The origin to be subtracted
  static void Transform(Polygon &polygon, const Vector2D &origin)
  {
    // Subtract origin from each vertex
    for (auto &p : polygon.Vertices)
      p -= origin;
  }

private:
  // Remove vertices from polygon keeping only vertices before given index
  static void RemoveVertices(Polygon &polygon, size_t end)
  {
    // Copy vertices to be kept to new vector
    std::vector<Vector2D> vertices(end);
    for (size_t i = 0; i < end; i++)
      vertices[i] = polygon.Vertices[i];

    // Overwrite vertices
    polygon.Vertices = vertices;
  }

  // Remove vertices from polygon (indices for removal assumed to be ordered)
  static void RemoveVertices(Polygon &polygon, const std::vector<size_t> remove)
  {
    // Copy vertices to be kept to new vector
    std::vector<Vector2D> vertices(polygon.Vertices.size() - remove.size());
    size_t k = 0;
    size_t l = 0;
    for (size_t i = 0; i < polygon.Vertices.size(); i++)
    {
      if (k < remove.size() && i == remove[k])
        k++;
      else
      {
        assert(l < vertices.size());
        vertices[l++] = polygon.Vertices[i];
      }
    }

    // Overwrite vertices
    polygon.Vertices = vertices;
  }
};

} // namespace DTCC

#endif
