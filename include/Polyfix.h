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
  static size_t MakeSimple(Polygon &polygon, double tol)
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

  /// Merge polygons. This creates a new polygon that covers the
  /// union of the two polygons and (as much as possible) respects
  /// the geometry of the two polygons.
  ///
  /// @param polygon0 First polygon
  /// @param polygon1 Second polygon
  /// @param tolerance Tolerance for connecting vertices and edges
  /// @return The merged polygon
  static Polygon
  Merge(const Polygon &polygon0, const Polygon &polygon1, double tol)
  {
    // Get number of vertices
    const size_t m = polygon0.Vertices.size();
    const size_t n = polygon1.Vertices.size();

    // Get all vertices
    std::vector<Point2D> vertices;
    vertices.reserve(m + n);
    for (const auto &p : polygon0.Vertices)
      vertices.push_back(p);
    for (const auto &p : polygon1.Vertices)
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

    // Find all pairwise connections between
    // edge i = (i0, i1) and edge j = (j0, j1)
    for (size_t i0 = 0; i0 < m; i0++)
    {
      const size_t i1 = edges[i0][0];
      for (size_t j0 = m; j0 < m + n; j0++)
      {
        const size_t j1 = edges[j0][0];

        // Find vertex-edge connections
        ConnectVertexEdge(i0, j0, j1, vertices, edges, tol);
        ConnectVertexEdge(i1, j0, j1, vertices, edges, tol);
        ConnectVertexEdge(j0, i0, i1, vertices, edges, tol);
        ConnectVertexEdge(j1, i0, i1, vertices, edges, tol);

        // Find edge-edge connections
        ConnectEdgeEdge(i0, i1, j0, j1, vertices, edges, tol);
      }
    }

    /*

    // Remove duplicate vertices
    numPoints = len(points)
    vertexMap = [i for i in range(numPoints)]
    removed = [false for i in range(numPoints)]
    for i in range(numPoints):
        for j in range(i+1, numPoints):
            if removed[i]: continue
            if DistancePointPoint(points[i], points[j]) < eps:
                #print('Merging:', i, j)
                edges[i] = edges[i] + edges[j]
                edges[j] = []
                vertexMap[j] = i
                removed[j] = true

    // Replace removed vertices in graph
    for i in range(numPoints):
        for j in range(len(edges[i])):
            edges[i][j] = vertexMap[edges[i][j]]

    // Remove duplicate edges in graph
    for i in range(len(edges)):
        newEdge = []
        for j in edges[i]:
            if j not in newEdge and i != j: newEdge.push_back(j)
        edges[i] = newEdge

    #for i, e in enumerate(edges):
    //    print(i, e)

    // Write point labels (and make sure they don't overlap)
    if plotting:
        H = 0.0075*max([Norm(p-q) for p in points for q in points])
        for i, p in enumerate(points):
            h = H
            for j, q in enumerate(points[:i]):
                if Norm(p - q) < eps:
                    h += 5*H
            text(p[0] + h, p[1] + H, str(i), va='bottom', ha='left')

    // Find first vertex by looking for an original edge that is to the
    // "right" of all points
    for i in range(m + n):

        // Skip if no outgoing edges
        if len(edges[i]) == 0: continue

        // Get the edge
        j = edges[i][0]
        u = points[j] - points[i]
        u /= Norm(u)

        // Check all points
        ok = true
        for k in range(numPoints):

            // Skip if removed
            if removed[k]: continue

            // Skip if on edge
            if k == i or k == j: continue

            // Check sin of angle (cross product)
            v = points[k] - points[i]
            v /= Norm(v)
            sin = u[0]*v[1] - u[1]*v[0]
            if sin < -eps:
                 ok = false
                break

        // Found first edge
        if ok:
            firstVertex = i
            nextVertex = j
            break

    // Keep track of visited vertices
    visited = [false for i in range(len(points))]
    visited[firstVertex] = true
    visited[nextVertex] = true

    // Initialize polygon
    vertices = [firstVertex, nextVertex]

    // Maximum number of step before failure
    maxNumSteps = 2*numPoints

    // Walk graph to build polygon counter-clockwise by picking
    // the right-most turn at each intersection
    for step in range(maxNumSteps):

        #print('Vertices:', vertices)

        // Get previous and current vertex
        i = len(vertices) - 1
        previousVertex = vertices[i - 1]
        currentVertex = vertices[i]

        // Get current edge(s)
        edge = edges[currentVertex]

        #print('Vertex %d: %s' % (currentVertex, str(edges[currentVertex])))

        // Find next vertex
        if len(edge) == 1:

            // If we only have one edge then follow it
            nextVertex = edge[0]

        else:

            // Get previous edge
            u = points[currentVertex] - points[previousVertex]
            d = Norm(u)
            u = u / d

            // Compute angles and distances for outgoing edges
            angles = []
            for k in edge:

                // Skip if already visited (if not first vertex)
                if k != firstVertex and visited[k]: continue

                // Skip if candidate edge intersects previous edges
                intersects = false
                e = (points[currentVertex], points[k])
                for l in range(1, i - 1):
                    f = (points[vertices[l]], points[vertices[l+1]])
                    if Intersects(e, f):
                        intersects = true
                        break
                if intersects: continue

                // Get new edge
                v = points[k] - points[currentVertex]
                d = Norm(v)
                v = v / d

                // Replace actual angle by cheaper but strictly increasing
                // function to avoid needing to call acos() or asin().
                sin = u[0]*v[1] - u[1]*v[0]
                cos = u[0]*v[0] + u[1]*v[1]
                if cos < -1.0 + eps: continue
                a = sin if cos >= 0.0 else (2.0-sin if sin > 0.0 else sin-2.0)
                angles.push_back((k, a, d))

            // If we have no more vertices to visit, take a step back
            if len(angles) == 0:
                print('No more vertices to visit, stepping back')
                del vertices[i]
                continue

            // Print angles
            #for angle in angles:
            //    print(angle)

            // Find smallest (right-most) angle. First priority is the angle
            // and second priority is the distance (pick closest). Note that
            // we add a small tolerance to ensure we get the closest vertex
            // if the vertices are on the same line.
            minAngle = angles[0]
            for angle in angles[1:]:
                if (angle[1] < minAngle[1] - eps) or \
                   (angle[1] < minAngle[1] + eps and angle[2] < minAngle[2]):
                    minAngle = angle

            // Pick next vertex
            nextVertex = minAngle[0]

        #print(currentVertex, "-->", nextVertex)
        #print('')

        // We are done if we return to the first vertex
        if nextVertex == firstVertex:
            #print('Back to first vertex')
            break

        // Add next vertex
        vertices.push_back(nextVertex)
        visited[nextVertex] = true

    // If merge failed, return convex hull
    if nextVertex != firstVertex:
        print('Merge failed, falling back to convex hull')
        points = [p for p in firstPolygon] + [p for p in secondPolygon]
        return ConvexHull(points)

    // Extract polygon points
    polygon = [points[i] for i in vertices]

    */

    return Polygon();
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

  // Connect vertex i to edge j = (j0, j1)
  static void ConnectVertexEdge(size_t i,
                                size_t j0,
                                size_t j1,
                                std::vector<Point2D> &vertices,
                                std::vector<std::vector<size_t>> &edges,
                                double tol)
  {
    // Avoid using sqrt for efficiency
    const double tol2 = tol * tol;

    // Get vertices
    const Point2D &p = vertices[i];
    const Point2D &q0 = vertices[j0];
    const Point2D &q1 = vertices[j1];

    // Connect vertices if close (create new edge)
    bool connected = false;
    if (Geometry::SquaredDistance2D(p, q0) < tol2)
    {
      edges[i].push_back(j0);
      edges[j0].push_back(i);
      connected = true;
    }
    if (Geometry::SquaredDistance2D(p, q1) < tol2)
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
    const double v2 = Geometry::SquaredNorm2D(v);
    if (v2 < Parameters::Epsilon)
      return;

    // Connect vertex to edge if close (project)
    if (Geometry::SquaredDistance2D(q0, q1, p) < tol2)
    {
      const Vector2D u(q0, p);
      const Point2D r = q0 + v * (Geometry::Dot2D(u, v) / v2);
      const size_t k = vertices.size();
      vertices.push_back(r);
      edges.push_back(std::vector<size_t>({i, j0, j1}));
      edges[i].push_back(k);
      edges[j0].push_back(k);
      edges[j1].push_back(k);
    }
  }

  // Connect edge (i0, i1) to edge (j0, j1)
  static void ConnectEdgeEdge(size_t i0,
                              size_t i1,
                              size_t j0,
                              size_t j1,
                              std::vector<Point2D> &vertices,
                              std::vector<std::vector<size_t>> &edges,
                              double tol)
  {
    // Avoid using sqrt for efficiency
    // const double tol2 = tol * tol;

    // Get vertices
    const Point2D &p0 = vertices[i0];
    const Point2D &p1 = vertices[i1];
    const Point2D &q0 = vertices[j0];
    const Point2D &q1 = vertices[j1];

    // Don't look for intersection if almost parallel
    const Vector2D u(p0, p1);
    const Vector2D v(q0, q1);
    if (std::abs(Geometry::Dot2D(u, v)) > 1.0 - Parameters::Epsilon)
      return;

    // Compute edge-edge intersection
    const Point2D r = Geometry::EdgeIntersection2D(p0, p1, q0, q1);

    // Connect edges to intersection if close
    if (Geometry::EdgeContains2D(p0, p1, r, tol) &&
        Geometry::EdgeContains2D(q0, q1, r, tol))
    {
      const size_t k = vertices.size();
      vertices.push_back(r);
      const int sp = Geometry::EdgeSign2D(p0, p1, r);
      const int sq = Geometry::EdgeSign2D(q0, q1, r);
      std::vector<size_t> kEdges{};
      if (sp == -1 || sp == 0)
      {
        edges[i0].push_back(k);
        kEdges.push_back(i0);
      }
      if (sp == 0 || sp == 1)
      {
        edges[i1].push_back(k);
        kEdges.push_back(i1);
      }
      if (sq == -1 || sq == 0)
      {
        edges[j0].push_back(k);
        kEdges.push_back(j0);
      }
      if (sq == 0 || sq == 1)
      {
        edges[j1].push_back(k);
        kEdges.push_back(j1);
      }
      edges.push_back(kEdges);
    }
  }
};

} // namespace DTCC

#endif
