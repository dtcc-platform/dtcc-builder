// Copyright (C) 2020-2021 Anders Logg, Anton J Olsson
// Licensed under the MIT License

#ifndef DTCC_POLYFIX_H
#define DTCC_POLYFIX_H

#include "Geometry.h"
#include "Logging.h"
#include "Timer.h"
#include "model/Polygon.h"

namespace DTCC_BUILDER
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
  static size_t make_closed(Polygon &polygon, double tol)
  {
    // Avoid using sqrt for efficiency
    const double tol2 = tol * tol;

    // Check each vertex against first vertex
    const Vector2D &p0 = Vector2D(polygon.vertices[0]);
    size_t end = 0;
    for (size_t i = 1; i < polygon.vertices.size() && end == 0; i++)
    {
      // Compute distance to first vertex
      const Vector2D &p = Vector2D(polygon.vertices[i]);
      const double d2 = Geometry::squared_distance_2d(p, p0);

      // Remove if distance is small
      if (i > 2 && d2 < tol2)
      {
        end = i;
        break;
      }
    }

    // Return if no vertices should be removed
    if (end == 0)
      return 0;

    // Remove vertices
    remove_vertices(polygon, end);

    return 1;
  }

  /// Make polygon counter-clockwise oriented.
  ///
  /// @param polygon The polygon
  /// @param tol Tolerance
  /// @return 0 if already counter-clockwise, 1 if modified
  static size_t make_oriented(Polygon &polygon)
  {
    // Return if already counter-clockwise
    if (Geometry::polygon_orientation_2d(polygon) == 0)
      return 0;

    // Reverse polygon
    std::reverse(polygon.vertices.begin(), polygon.vertices.end());

    return 1;
  }

  /// Merge polygon vertices.
  ///
  /// @param polygon The polygon
  /// @param tol Tolerance for small distance
  /// @return 0 if already simple, 1 if modified
  static size_t merge_vertices(Polygon &polygon, double tol)
  {
    // Avoid using sqrt for efficiency
    const double tol2 = tol * tol;

    // vertices to be removed
    std::vector<size_t> remove{};

    // Check each edge
    const size_t num_vertices = polygon.vertices.size();
    size_t i0 = 0;
    for (size_t i = 1; i < num_vertices; i++)
    {
      // Get current vertex
      const Vector2D &p0 = Vector2D(polygon.vertices[i0]);
      const Vector2D &p = Vector2D(polygon.vertices[i]);

      // Compute distance
      const double d2 = Geometry::squared_distance_2d(p0, p);

      // Remove if distance is small
      if (d2 < tol2)
        remove.push_back(i);
      else
        i0 = i;
    }

    // Return if no vertices should be removed
    if (remove.empty())
      return 0;

    // Remove vertices
    remove_vertices(polygon, remove);

    return 1;
  }

  /// Merge polygon edges.
  ///
  /// @param polygon The polygon
  /// @param tol Tolerance for small angle (sin of angle)
  /// @return 0 if already simple, 1 if modified
  static size_t merge_edges(Polygon &polygon, double tol)
  {
    // Avoid using sqrt for efficiency
    const double tol2 = tol * tol;

    // vertices to be removed
    std::vector<size_t> remove{};

    // Check each edge
    const size_t num_vertices = polygon.vertices.size();
    for (size_t i = 0; i < num_vertices; i++)
    {
      // Get previous, current and next points
      const Vector2D &p0 =
          Vector2D(polygon.vertices[(i + num_vertices - 1) % num_vertices]);
      const Vector2D &p1 = Vector2D(polygon.vertices[i]);
      const Vector2D &p2 = Vector2D(polygon.vertices[(i + 1) % num_vertices]);

      // Compute edges and dot products
      const Vector2D u = p1 - p0;
      const Vector2D v = p2 - p1;
      const double u2 = Geometry::dot_2d(u, u);
      const double v2 = Geometry::dot_2d(v, v);
      const double uv = Geometry::dot_2d(u, v);

      // Remove if angle is small
      if (uv * uv > (1.0 - tol2) * u2 * v2)
        remove.push_back(i);
    }

    // Return if no vertices should be removed
    if (remove.empty())
      return 0;

    // Remove vertices
    remove_vertices(polygon, remove);

    return 1;
  }

  /// transform polygon by subtracting given origin.
  ///
  /// @param polygon The polygon
  /// @param origin The origin to be subtracted
  static void transform(Polygon &polygon, const Vector2D &origin)
  {
    transform(polygon.vertices, origin);
  }

  /// transform vertices by subtracting given origin.
  ///
  /// @param polygon The polygon
  /// @param origin The origin to be subtracted
  static void transform(std::vector<Vector2D> &vertices, const Vector2D &origin)
  {
    // Subtract origin from each vertex
    Vector2D _origin(origin);
    for (auto &p : vertices)
      p -= _origin;
  }

  /// Merge polygons. This creates a new polygon that covers the
  /// union of the two polygons and (as much as possible) respects
  /// the geometry of the two polygons.
  ///
  /// @param polygon0 First polygon
  /// @param polygon1 Second polygon
  /// @param tolerance Tolerance for connecting vertices and edges
  /// @return The merged polygon
  static Polygon
  merge_polygons(const Polygon &polygon0, const Polygon &polygon1, double tol)
  {
    Timer timer("merge_polygons");

    // Avoid using sqrt for efficiency
    // const double tol2 = tol * tol;
    const double eps = Constants::epsilon;
    const double eps2 = eps * eps;
    const double tol2 = tol * tol;

    // Get number of vertices
    const size_t m = polygon0.vertices.size();
    const size_t n = polygon1.vertices.size();

    // Get all vertices
    std::vector<Vector2D> vertices;
    vertices.reserve(m + n);
    for (const auto &p : polygon0.vertices)
      vertices.push_back(p);
    for (const auto &p : polygon1.vertices)
      vertices.push_back(p);

    // Create directed graph of edges
    std::vector<std::vector<size_t>> edges;
    edges.reserve(m + n);
    for (size_t i = 0; i < m; i++)
    {
      const std::vector<size_t> edge = {(i + 1) % m};
      edges.push_back(edge);
    }
    for (size_t i = 0; i < n; i++)
    {
      const std::vector<size_t> edge = {(i + 1) % n + m};
      edges.push_back(edge);
    }

    // print_edges(edges);

    // find all pairwise connections between
    // edge i = (i0, i1) and edge j = (j0, j1)
    for (size_t i0 = 0; i0 < m; i0++)
    {
      const size_t i1 = edges[i0][0];
      for (size_t j0 = m; j0 < m + n; j0++)
      {
        const size_t j1 = edges[j0][0];

        // find vertex-edge connections
        connect_vertex_edge(i0, j0, j1, vertices, edges, tol);
        connect_vertex_edge(i1, j0, j1, vertices, edges, tol);
        connect_vertex_edge(j0, i0, i1, vertices, edges, tol);
        connect_vertex_edge(j1, i0, i1, vertices, edges, tol);

        // find edge-edge connections
        connect_edge_edge(i0, i1, j0, j1, vertices, edges, tol);
      }
    }

    // print_edges(edges);

    // Remove duplicate vertices
    assert(vertices.size() == edges.size());
    const size_t num_vertices = vertices.size();
    std::vector<size_t> vertex_map(num_vertices);
    std::vector<bool> removed(num_vertices);
    for (size_t i = 0; i < num_vertices; i++)
    {
      vertex_map[i] = i;
      removed[i] = false;
    }
    for (size_t i = 0; i < num_vertices; i++)
    {
      for (size_t j = i + 1; j < num_vertices; j++)
      {
        if (removed[i])
          continue;
        if (Geometry::squared_distance_2d(vertices[i], vertices[j]) < tol2)
        {
          for (const auto k : edges[j])
            edges[i].push_back(k);
          edges[j].clear();
          vertex_map[j] = i;
          removed[j] = true;
        }
      }
    }

    // Replace removed vertices in graph
    for (size_t i = 0; i < edges.size(); i++)
      for (auto &edge : edges[i])
        edge = vertex_map[edge];

    // Remove duplicate edges in graph
    for (size_t i = 0; i < edges.size(); i++)
    {
      std::vector<size_t> new_edges{};
      for (const auto &edge : edges[i])
      {
        if (edge != i && std::find(new_edges.begin(), new_edges.end(), edge) ==
                             new_edges.end())
          new_edges.push_back(edge);
      }
      edges[i] = new_edges;
    }

    // print_edges(edges);

    // find first vertex by looking for an original edge that is to the
    // "right" of all points
    size_t first_vertex{}, next_vertex{};
    for (size_t i = 0; i < m + n; i++)
    {
      // Skip if no outgoing edges
      if (edges[i].empty())
        continue;

      // Get the edge
      const size_t j = edges[i][0];
      const Vector2D u(vertices[i], vertices[j]);
      const double u2 = Geometry::squared_norm_2d(u);

      // Check all points
      bool ok = true;
      for (size_t k = 0; k < num_vertices; k++)
      {
        // Skip if removed
        if (removed[k])
          continue;

        // Skip if on edge
        if (k == i || k == j)
          continue;

        // Check sin of angle (cross product)
        const Vector2D v(vertices[i], vertices[k]);
        const double v2 = Geometry::squared_norm_2d(v);
        const double sin = u.x * v.y - u.y * v.x;
        if (sin < 0.0 && sin * sin > eps2 * u2 * v2)
        {
          ok = false;
          break;
        }
      }

      // Found first edge
      if (ok)
      {
        first_vertex = i;
        next_vertex = j;
        break;
      }
    }

    // Keep track of visited vertices
    std::vector<bool> visited(num_vertices);
    std::fill(visited.begin(), visited.end(), false);
    visited[first_vertex] = true;
    visited[next_vertex] = true;

    // Initialize polygon
    std::vector<size_t> polygon;
    polygon.push_back(first_vertex);
    polygon.push_back(next_vertex);

    // Maximum number of step before failure
    const size_t max_num_steps = 2 * num_vertices;

    // Walk graph to build polygon counter-clockwise by picking
    // the right-most turn at each intersection
    for (size_t step = 0; step < max_num_steps; step++)
    {
      // Get previous and current vertex
      const size_t i = polygon.size() - 1;
      const size_t previous_vertex = polygon[i - 1];
      const size_t current_vertex = polygon[i];

      // Get current edge(s)
      const std::vector<size_t> &edge = edges[current_vertex];

      // find next vertex
      assert(!edge.empty());
      if (edge.size() == 1)
      {
        // If we only have one edge then follow it
        next_vertex = edge[0];
      }
      else
      {
        // Get previous edge
        const Vector2D u(vertices[previous_vertex], vertices[current_vertex]);

        // Compute angles and distances for outgoing edges
        std::vector<std::tuple<size_t, double, double>> candidates{};
        for (const auto k : edge)
        {
          // Skip if already visited (if not first vertex)
          if (k != first_vertex && visited[k])
            continue;

          // Skip if candidate edge intersects a previous edge
          bool intersects = false;
          const auto e = std::make_pair(vertices[current_vertex], vertices[k]);
          for (size_t l = 1; l < i - 1; l++)
          {
            const auto f =
                std::make_pair(vertices[polygon[l]], vertices[polygon[l + 1]]);
            if (Geometry::intersects_2d(e.first, e.second, f.first, f.second))
            {
              intersects = true;
              break;
            }
          }
          if (intersects)
            continue;

          // Get new edge
          const Vector2D v(vertices[current_vertex], vertices[k]);
          const double v2 = Geometry::squared_norm_2d(v);

          // Compute vector angle and add candidate edge
          const double a = Geometry::vector_angle_2d(u, v);
          if (std::abs(a) > 2.0 - eps)
            continue;
          candidates.push_back(std::make_tuple(k, a, v2));
          // std::cout << "Adding k = " << k << " a = " << a << std::endl;
        }

        // If we have no more vertices to visit, take a step back
        if (candidates.empty())
        {
          polygon.pop_back();
          continue;
        }

        // find smallest (right-most) angle. First priority is the angle
        // and second priority is the distance (pick closest). Note that
        // we add a small tolerance to ensure we get the closest vertex
        // if the vertices are on the same line.
        auto min_candidate = candidates[0];
        for (size_t k = 1; k < candidates.size(); k++)
        {
          const double angle = std::get<1>(candidates[k]);
          const double distance = std::get<2>(candidates[k]);
          const double min_angle = std::get<1>(min_candidate);
          const double min_distance = std::get<2>(min_candidate);
          if ((angle < min_angle - 0.01) ||
              (angle < min_angle + 0.01 && distance < min_distance))
          {
            min_candidate = candidates[k];
          }
        }

        // Pick next vertex
        next_vertex = std::get<0>(min_candidate);
      }

      // We are done if we return to the first vertex
      if (next_vertex == first_vertex)
        break;

      // Add next vertex
      polygon.push_back(next_vertex);
      visited[next_vertex] = true;
    }

    // If merge failed, return convex hull
    if (next_vertex != first_vertex)
    {
      warning("Polygon merge failed, falling back to convex hull.");
      vertices.clear();
      for (const auto &p : polygon0.vertices)
        vertices.push_back(p);
      for (const auto &p : polygon1.vertices)
        vertices.push_back(p);
      return Geometry::convex_hull_2d(vertices);
    }

    // print_polygon(polygon);

    // Create polygon
    Polygon _polygon{};
    _polygon.vertices.reserve(polygon.size());
    for (size_t i = 0; i < polygon.size(); i++)
      _polygon.vertices.push_back(vertices[polygon[i]]);

    return _polygon;
  }

private:
  // Remove vertices from polygon keeping only vertices before given index
  static void remove_vertices(Polygon &polygon, size_t end)
  {
    // Copy vertices to be kept to new vector
    std::vector<Vector2D> vertices(end);
    for (size_t i = 0; i < end; i++)
      vertices[i] = polygon.vertices[i];

    // Overwrite vertices
    polygon.vertices = vertices;
  }

  // Remove vertices from polygon (indices for removal assumed to be ordered)
  static void remove_vertices(Polygon &polygon,
                              const std::vector<size_t> &remove)
  {
    // Copy vertices to be kept to new vector
    std::vector<Vector2D> vertices(polygon.vertices.size() - remove.size());
    size_t k = 0;
    size_t l = 0;
    for (size_t i = 0; i < polygon.vertices.size(); i++)
    {
      if (k < remove.size() && i == remove[k])
        k++;
      else
      {
        assert(l < vertices.size());
        vertices[l++] = polygon.vertices[i];
      }
    }

    // Overwrite vertices
    polygon.vertices = vertices;
  }

  // Connect vertex i to edge j = (j0, j1)
  static void connect_vertex_edge(size_t i,
                                  size_t j0,
                                  size_t j1,
                                  std::vector<Vector2D> &vertices,
                                  std::vector<std::vector<size_t>> &edges,
                                  double tol)
  {
    // Avoid using sqrt for efficiency
    const double tol2 = tol * tol;
    const double eps2 = Constants::epsilon * Constants::epsilon;

    // Get vertices
    const Vector2D &p = vertices[i];
    const Vector2D &q0 = vertices[j0];
    const Vector2D &q1 = vertices[j1];

    // Connect vertices if close (create new edge)
    bool connected = false;
    if (Geometry::squared_distance_2d(p, q0) < tol2)
    {
      edges[i].push_back(j0);
      edges[j0].push_back(i);
      connected = true;
    }
    if (Geometry::squared_distance_2d(p, q1) < tol2)
    {
      edges[i].push_back(j1);
      edges[j1].push_back(i);
      connected = true;
    }

    // Don't connect vertex to edge if already connected
    if (connected)
      return;

    // Don't connect vertex to edge if zero length
    Vector2D v(q0, q1);
    const double v2 = Geometry::squared_norm_2d(v);
    if (v2 < eps2)
      return;

    // Connect vertex to edge if close (project)
    if (Geometry::squared_distance_2d(q0, q1, p) < tol2)
    {
      const Vector2D u(q0, p);
      const Vector2D r = q0 + v * (Geometry::dot_2d(u, v) / v2);
      const size_t k = vertices.size();
      vertices.push_back(r);
      edges.push_back(std::vector<size_t>({i, j0, j1}));
      edges[i].push_back(k);
      edges[j0].push_back(k);
      edges[j1].push_back(k);
    }
  }

  // Connect edge (i0, i1) to edge (j0, j1)
  static void connect_edge_edge(size_t i0,
                                size_t i1,
                                size_t j0,
                                size_t j1,
                                std::vector<Vector2D> &vertices,
                                std::vector<std::vector<size_t>> &edges,
                                double tol)
  {
    // Get vertices (don't use references since vector might increase in size)
    assert(i0 < vertices.size());
    assert(i1 < vertices.size());
    assert(j0 < vertices.size());
    assert(j1 < vertices.size());
    const Vector2D p0 = vertices[i0];
    const Vector2D p1 = vertices[i1];
    const Vector2D q0 = vertices[j0];
    const Vector2D q1 = vertices[j1];

    // Don't look for intersection if almost parallel
    const Vector2D u(p0, p1);
    const Vector2D v(q0, q1);
    const double uv = Geometry::dot_2d(u, v);
    const double u2 = Geometry::squared_norm_2d(u);
    const double v2 = Geometry::squared_norm_2d(v);
    const double eps = Constants::epsilon;
    if (uv * uv > (1.0 - eps) * (1.0 - eps) * u2 * v2)
      return;

    // Compute edge-edge intersection
    const Vector2D r = Geometry::edge_intersection_2d(p0, p1, q0, q1);

    // Connect edges to intersection if close
    if (Geometry::edge_contains_2d(p0, p1, r, tol) &&
        Geometry::edge_contains_2d(q0, q1, r, tol))
    {
      const size_t k = vertices.size();
      vertices.push_back(r);
      const int sp = Geometry::edge_sign_2d(p0, p1, r);
      const int sq = Geometry::edge_sign_2d(q0, q1, r);
      std::vector<size_t> k_edges{};
      if (sp == -1 || sp == 0)
      {
        edges[i0].push_back(k);
        k_edges.push_back(i0);
      }
      if (sp == 0 || sp == 1)
      {
        edges[i1].push_back(k);
        k_edges.push_back(i1);
      }
      if (sq == -1 || sq == 0)
      {
        edges[j0].push_back(k);
        k_edges.push_back(j0);
      }
      if (sq == 0 || sq == 1)
      {
        edges[j1].push_back(k);
        k_edges.push_back(j1);
      }
      edges.push_back(k_edges);
    }
  }

  // Used for debugging
  static void print_edges(const std::vector<std::vector<size_t>> &edges)
  {
    std::stringstream string_builder{};
    string_builder << "Edges: [";
    for (size_t i = 0; i < edges.size(); i++)
    {
      if (i > 0)
        string_builder << ", ";
      string_builder << "[";
      for (size_t j = 0; j < edges[i].size(); j++)
      {
        if (j > 0)
          string_builder << ", ";
        string_builder << edges[i][j];
      }
      string_builder << "]";
    }
    string_builder << "]" << std::endl;
    info(string_builder.str());
  }

  // Used for debugging
  static void print_polygon(const std::vector<size_t> &polygon)
  {
    std::stringstream string_builder{};
    string_builder << "Polygon: [";
    for (size_t i = 0; i < polygon.size(); i++)
    {
      if (i > 0)
        string_builder << ", ";
      string_builder << polygon[i];
    }
    string_builder << "]" << std::endl;
    info(string_builder.str());
  }
};

} // namespace DTCC_BUILDER

#endif
