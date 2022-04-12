// Copyright (C) 2022 Anders Logg
// Licensed under the MIT License

#ifndef GEOS_H
#define GEOS_H

#include <geos_c.h>

#include "Logging.h"
#include "Polygon.h"
#include "Timer.h"

// GEOS message handler
static void geos_msg_handler(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  vprintf(fmt, ap);
  va_end(ap);
}

namespace DTCC
{
/// This class provides algorithms for working with polygons based on
/// the GEOS library (libgeos).
///
///
/// Note: GEOS::Init() needs to be called before using any of the
/// methods and GEOS::Finish() needs to be called to clean up
/// resources.
class GEOS
{
public:
  /// Initialize GEOS, needs to be called before any other methods
  static void Init()
  {
    Info("Initializing GEOS");
    initGEOS(geos_msg_handler, geos_msg_handler);
  }

  /// Finish GEOS, needs to be called to clean up resources
  static void Finish()
  {
    Info("Cleaning up GEOS");
    finishGEOS();
  }

  /// Print GEOS geometry to terminal (in WKT format)
  static void Print(GEOSGeometry *geometry)
  {
    // Print geometry as WKT
    GEOSWKTWriter *writer = GEOSWKTWriter_create();
    GEOSWKTWriter_setTrim(writer, 1);
    char *out = GEOSWKTWriter_write(writer, geometry);
    Info(out);

    // Cleanup
    GEOSWKTWriter_destroy(writer);
    GEOSFree(out);
  }

  /// Create GEOS geoemetry from Polygon
  static GEOSGeometry *CreateGeometry(const Polygon &polygon)
  {
    // Check that polygon is not empty
    const size_t n = polygon.Vertices.size();
    if (n == 0)
      Error("Unable to create GEOS polygon; polygon has no vertices");

    // Create coordinate sequence and set coordinates. Note that first
    // vertex needs to be included both at beginning and end.
    GEOSCoordSequence *sequence = GEOSCoordSeq_create(n + 1, 2);
    for (size_t i = 0; i <= n; i++)
    {
      GEOSCoordSeq_setX(sequence, i, polygon.Vertices[i % n].x);
      GEOSCoordSeq_setY(sequence, i, polygon.Vertices[i % n].y);
    }

    // Create polygon from coordinate sequence
    GEOSGeometry *ring = GEOSGeom_createLinearRing(sequence);
    GEOSGeometry *geometry = GEOSGeom_createPolygon(ring, 0, 0);

    return geometry;
  }

  /// Create Polygon from GEOS geometry
  static Polygon CreatePolygon(const GEOSGeometry *geometry)
  {
    Polygon p;
    return p;
  }

  /// Merge polygons. This creates a new polygon that covers the
  /// union of the two polygons and (as much as possible) respects
  /// the geometry of the two polygons.
  ///
  /// This version using libgeos to computed the merged polygon.
  ///
  /// @param polygon0 First polygon
  /// @param polygon1 Second polygon
  /// @param tolerance Tolerance for connecting vertices and edges
  /// @return The merged polygon
  static Polygon
  MergePolygons(const Polygon &polygon0, const Polygon &polygon1, double tol)
  {
    Timer timer("MergePolygons (GEOS)");

    // FIXME: Shouldn't be here
    const double EPS = 0.01;

    // Create GEOS geometries
    GEOSGeometry *A = CreateGeometry(polygon0);
    GEOSGeometry *B = CreateGeometry(polygon1);

    // Simplify geometries
    GEOSGeometry *_A = GEOSSimplify(A, tol);
    GEOSGeometry *_B = GEOSSimplify(B, tol);

    // Compute union
    GEOSGeometry *C = GEOSUnionPrec(A, B, EPS);
    GEOSGeometry *_C = GEOSSimplify(C, tol);

    // Geometry to be created
    GEOSGeometry *geometry{};

    // Accept if polygon (not multi-polygon) and good quality
    if (IsValid(_C, tol))
    {
      Info("Using union");
      geometry = _C;
    }

    // If not acceptable, try gradually increase tolerance
    if (geometry == 0)
    {
      double _tol = tol;
      const size_t maxiter = 3;
      for (size_t k = 0; k < maxiter; k++)
      {
        // Compute vertex projections
        std::vector<Vector2D> projections{};
        ComputeVertexProjections(_A, _B, tol, projections);
        ComputeVertexProjections(_B, _A, tol, projections);

        // Check that we have at least three vertices
        if (projections.size() >= 3)
        {
          // Create GEOS geometry for projections
          GEOSCoordSequence *sequence =
              GEOSCoordSeq_create(projections.size(), 2);
          for (size_t i = 0; i < projections.size(); i++)
          {
            GEOSCoordSeq_setX(sequence, i, projections[i].x);
            GEOSCoordSeq_setY(sequence, i, projections[i].y);
          }
          GEOSGeometry *string = GEOSGeom_createLineString(sequence);

          // Compute convex hull of projections
          GEOSGeometry *P = GEOSConvexHull(string);

          // Compute union
          GEOSGeometry *AB = GEOSUnionPrec(_A, _B, EPS);
          GEOSGeometry *D = GEOSUnionPrec(AB, P, EPS);
          GEOSGeometry *_D = GEOSSimplify(D, tol);

          // Accept if single geonetry and good quality
          if (IsValid(_D, tol))
          {
            Info("Using extended union");
            geometry = _D;
            break;
          }
        }

        // Increase tolerance
        _tol *= 2;
        std::cout << "Increasing tolerance: tol = " << _tol << std::endl;
      }
    }

    // If not acceptable, try merging convex hulls
    if (geometry == 0)
    {
      GEOSGeometry *HA = GEOSConvexHull(_A);
      GEOSGeometry *HB = GEOSConvexHull(_B);
      GEOSGeometry *D = GEOSUnionPrec(HA, HB, EPS);
      GEOSGeometry *_D = GEOSSimplify(D, tol);

      // Accept if single geonetry and good quality
      if (IsValid(_D, tol))
      {
        Info("Using union of convex hulls");
        geometry = _D;
      }
    }

    // If not acceptable, use convex hull
    if (geometry == 0)
    {
      GEOSGeometry *_AB = GEOSUnionPrec(_A, _B, EPS);
      GEOSGeometry *H = GEOSConvexHull(_AB);
      GEOSGeometry *_H = GEOSSimplify(H, tol);
      geometry = _H;
    }

    // Convert to polygon
    Print(geometry);
    Polygon polygon = CreatePolygon(geometry);

    // Cleanup
    GEOSGeom_destroy(A);
    GEOSGeom_destroy(B);
    GEOSGeom_destroy(C);
    GEOSGeom_destroy(_A);
    GEOSGeom_destroy(_B);
    GEOSGeom_destroy(_C);

    return polygon;
  }

private:
  // Check if geometry is valid
  static bool IsValid(const GEOSGeometry *geometry, double tol)
  {
    // Get minimum clearance
    double q{};
    GEOSMinimumClearance(geometry, &q);

    // Get type
    const bool isPolygon = GEOSGeomTypeId(geometry) == GEOS_POLYGON;

    return isPolygon && q > tol;
  }

  // Project vertices of polygon A on edges of polygon B if close
  static void ComputeVertexProjections(const GEOSGeometry *A,
                                       const GEOSGeometry *B,
                                       double tol,
                                       std::vector<Vector2D> &projections)
  {
    std::cout << "Computing vertex projections" << std::endl;

    // Get vertices
    const GEOSGeometry *rA = GEOSGetExteriorRing(A);
    const GEOSGeometry *rB = GEOSGetExteriorRing(B);
    const GEOSCoordSequence *vA = GEOSGeom_getCoordSeq(rA);
    const GEOSCoordSequence *vB = GEOSGeom_getCoordSeq(rB);
    const int nA = GEOSGetNumCoordinates(A);
    const int nB = GEOSGetNumCoordinates(B);

    // Iterate over vertices in A
    for (int i = 0; i < (nA - 1); i++)
    {
      // Get vertex
      Vector2D q{};
      GEOSCoordSeq_getXY(vA, i, &q.x, &q.y);

      // Iterate over edges in B
      for (int j = 0; j < (nB - 1); j++)
      {
        // Get edge
        Vector2D p0{}, p1{};
        GEOSCoordSeq_getXY(vB, j, &p0.x, &p0.y);
        GEOSCoordSeq_getXY(vB, j + 1, &p1.x, &p1.y);

        // Project point to edge
        Vector2D u = q - p0;
        Vector2D v = p1 - p0;
        Vector2D p = p0 + v * Geometry::Dot2D(u, v) / Geometry::Dot2D(v, v);

        // Check whether projected point is inside segment. Check either
        // x or y coordinates depending on which is largest (most stable)
        bool inside{};
        if (std::abs(v.x) > std::abs(v.y))
          inside = std::min(p0.x, p1.x) <= p.x && p.x <= std::max(p0.x, p1.x);
        else
          inside = std::min(p0.y, p1.y) <= p.y && p.y <= std::max(p0.y, p1.y);

        // Use either project or closest vertex
        if (inside)
        {
          const double d = sqrt(Geometry::Dot2D(q - p, q - p));
          if (d < tol)
          {
            projections.push_back(p);
            projections.push_back(q);
          }
          else
          {
            const double d0 = sqrt(Geometry::Dot2D(q - p0, q - p0));
            const double d1 = sqrt(Geometry::Dot2D(q - p1, q - p1));
            if (d0 < d1 && d0 < tol)
            {
              projections.push_back(p0);
              projections.push_back(q);
            }
            else if (d1 < tol)
            {
              projections.push_back(p1);
              projections.push_back(q);
            }
          }
        }
      }
    }
  }

  /*
    # Simplify geometries
    A = simplify(A, TOL)
    B = simplify(B, TOL)

    # Compute union
    C = union(A, B, grid_size=EPS)
    C = simplify(C, TOL)

    # Accept if single geometry and good quality
    if IsValid(C):
    print('Using union')
    return C

    # Gradually increase tolerance for mergin
    tol = TOL
    maxiter = 3
    for k in range(maxiter):

    # Compute vertex projections
    pAB = ComputeVertexProjections(A, B, tol)
    pBA = ComputeVertexProjections(B, A, tol)
    projections = pAB + pBA

    # Check that we have at least vertices
    if len(projections) >= 3:

    # Compute convex hull of projections
    P = convex_hull(polygons(projections))

    # Compute union
    C = union_all([A, B, P], grid_size=EPS)
    C = simplify(C, TOL)

    # Accept if single geonetry and good quality
    if IsValid(C):
    print('Using extended union')
    return C

    # Increase tolerance
    tol *= 2
    print('Increasing tolerance: tol =', tol)

    # Try merging convex hulls
    A = convex_hull(A)
    B = convex_hull(B)
    C = union(A, B, grid_size=EPS)
    C = simplify(C, TOL)
    if IsValid(C):
    print('Using union of convex hulls')
    return C

    # Return convex hull
    C = union(A, B)
    C = convex_hull(C)
    C = simplify(C, TOL)

    print('Using convex hull')
    return C
  */
};

} // namespace DTCC

#endif
