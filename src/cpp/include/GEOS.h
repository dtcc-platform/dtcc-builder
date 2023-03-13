// Copyright (C) 2022 Anders Logg
// Licensed under the MIT License

#ifndef GEOS_H
#define GEOS_H

#include <geos_c.h>
#include <stdarg.h>

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

namespace DTCC_BUILDER
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
    info("GEOS: Initializing");
    initGEOS(geos_msg_handler, geos_msg_handler);
  }

  /// Finish GEOS, needs to be called to clean up resources
  static void Finish()
  {
    info("GEOS: Cleaning up");
    finishGEOS();
  }

  /// Print GEOS geometry to terminal (in WKT format)
  static void Print(GEOSGeometry *geometry)
  {
    // Print geometry as WKT
    GEOSWKTWriter *writer = GEOSWKTWriter_create();
    GEOSWKTWriter_setTrim(writer, 1);
    char *out = GEOSWKTWriter_write(writer, geometry);
    info(out);

    // Cleanup
    GEOSWKTWriter_destroy(writer);
    GEOSFree(out);
  }

  /// Create GEOS geoemetry from Polygon. This creates a new geometry
  /// and caller is responsible for calling GEOSGeom_destroy().
  static GEOSGeometry *CreateGeometry(const Polygon &polygon)
  {
    // Check that polygon is not empty
    const size_t n = polygon.Vertices.size();
    if (n == 0)
      error("Unable to create GEOS polygon; polygon has no vertices");

    // Create coordinate sequence and set coordinates. Note that first
    // vertex needs to be included both at beginning and end.
    GEOSCoordSequence *sequence = GEOSCoordSeq_create(n + 1, 2);
    for (size_t i = 0; i <= n; i++)
    {
      GEOSCoordSeq_setX(sequence, i, polygon.Vertices[i % n].x);
      GEOSCoordSeq_setY(sequence, i, polygon.Vertices[i % n].y);
    }

    // Create polygon from coordinate sequence. Note that ownership
    // is transferred to the new objects.
    GEOSGeometry *ring = GEOSGeom_createLinearRing(sequence);
    GEOSGeometry *geometry = GEOSGeom_createPolygon(ring, 0, 0);

    return geometry;
  }

  /// Create Polygon from GEOS geometry
  static Polygon CreatePolygon(const GEOSGeometry *geometry)
  {
    // Get coordinate sequence
    const GEOSGeometry *ring = GEOSGetExteriorRing(geometry);
    const GEOSCoordSequence *sequence = GEOSGeom_getCoordSeq(ring);
    const int n = GEOSGetNumCoordinates(ring);

    // Create polygon. Note that we skip the last (duplicate) vertex.
    Polygon polygon;
    for (int i = 0; i < (n - 1); i++)
    {
      Point2D p{};
      GEOSCoordSeq_getXY(sequence, i, &p.x, &p.y);
      polygon.Vertices.push_back(p);
    }

    return polygon;
  }

  /// Simplify a polygon using the Douglas-Peucker algorithm. This
  /// function is a wrapper around the function GEOSSimplify(), but
  /// it guarantees that the result is a Polygon. If it fails to simplify
  /// the polygon, it will return the original Polygon.
  ///
  /// @param geometry The geometry to simplify
  /// @param tol The tolerance for the simplification
  /// @param max_iter The maximum number of iterations to try to simplify
  /// @return The simplified polygon
  static GEOSGeometry *
  SimplifyPolygon(const GEOSGeometry *geometry, double tol, size_t max_iter = 3)
  {

    for (size_t i = 0; i < max_iter; i++)
    {
      GEOSGeometry *simplified = GEOSSimplify(geometry, tol);
      if (GEOSGeomTypeId(simplified) == GEOS_POLYGON)
      {
        return simplified;
      }
      else
      {
        GEOSGeom_destroy(simplified);
        tol /= 2;
      }
    }
    GEOSGeometry *simplified = GEOSGeom_clone(geometry);

    warning("Unable to Simplify polygon; revert to originial");
    return simplified;
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

    // Set geometric precision
    const double EPS = Parameters::Epsilon;

    // Create A, simplify, and cleanup
    GEOSGeometry *_A = CreateGeometry(polygon0);
    GEOSGeometry *A = SimplifyPolygon(_A, tol);
    GEOSGeom_destroy(_A);
    _A = 0;

    // Create B, simplify, and cleanup
    GEOSGeometry *_B = CreateGeometry(polygon1);
    GEOSGeometry *B = SimplifyPolygon(_B, tol);
    GEOSGeom_destroy(_B);
    _B = 0;

    // Compute union C, simplify, and cleanup
    GEOSGeometry *_C = GEOSUnionPrec(A, B, EPS);
    GEOSGeometry *C = SimplifyPolygon(_C, tol);
    GEOSGeom_destroy(_C);
    _C = 0;

    // A note regarding memory allocation and cleanup... At this point
    // we have allocated three geometries A, B, C that need to be
    // cleaned up below. Any intermediate geometries allocated below
    // should be cleaned up close to the point of allocation.

    // Accept if valid
    if (IsValid(C, tol))
    {
      // info("Using union");
    }
    else
    {
      GEOSGeom_destroy(C);
      C = 0;
    }

    // Print(A);
    // Print(B);

    // If not accepted, increase tolerance
    if (C == 0)
    {
      double _tol = tol;
      const size_t maxiter = 3;
      for (size_t k = 0; k < maxiter; k++)
      {
        // Compute vertex projections
        std::vector<Vector2D> projections{};
        ComputeVertexProjections(A, B, _tol, projections);
        ComputeVertexProjections(B, A, _tol, projections);

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

          // Cleanup (note that string takes over sequence)
          GEOSGeom_destroy(string);

          // Compute union
          GEOSGeometry *AB = GEOSUnionPrec(A, B, EPS);
          GEOSGeometry *_C = GEOSUnionPrec(AB, P, EPS);
          C = SimplifyPolygon(_C, tol);

          // Cleanup
          GEOSGeom_destroy(P);
          GEOSGeom_destroy(AB);
          GEOSGeom_destroy(_C);

          // Accept if valid
          if (IsValid(C, tol))
          {
            // info("Using extended union");
            break;
          }
          else
          {
            GEOSGeom_destroy(C);
            C = 0;
          }
        }

        // Increase tolerance
        _tol *= 2;
      }
    }

    // If not acceptable, try union of convex hulls
    if (C == 0)
    {
      // Compute union of convex hulls
      GEOSGeometry *HA = GEOSConvexHull(A);
      GEOSGeometry *HB = GEOSConvexHull(B);
      GEOSGeometry *_C = GEOSUnionPrec(HA, HB, EPS);
      C = SimplifyPolygon(_C, tol);

      // Cleanup
      GEOSGeom_destroy(HA);
      GEOSGeom_destroy(HB);
      GEOSGeom_destroy(_C);

      // Accept if valid
      if (IsValid(C, tol))
      {
        warning("GEOS: Falling back to union of convex hulls");
      }
      else
      {
        GEOSGeom_destroy(C);
        C = 0;
      }
    }

    // If not acceptable, use convex hull of union
    if (C == 0)
    {
      warning("GEOS: Falling back to convex hull of union");

      // Compute convex hull of union
      GEOSGeometry *AB = GEOSUnionPrec(A, B, EPS);
      GEOSGeometry *_C = GEOSConvexHull(AB);
      C = SimplifyPolygon(_C, tol);

      // Cleanup
      GEOSGeom_destroy(AB);
      GEOSGeom_destroy(_C);
    }

    // Convert to polygon
    assert(C);
    // Print(C);
    Polygon polygon = CreatePolygon(C);

    // Cleanup
    assert(A);
    assert(B);
    assert(C);
    GEOSGeom_destroy(A);
    GEOSGeom_destroy(B);
    GEOSGeom_destroy(C);

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
    // Get vertices
    if (GEOSGeomTypeId(A) != GEOS_POLYGON)
      info("A is not a polygon");
    if (GEOSGeomTypeId(B) != GEOS_POLYGON)
    {
      info("B is not a polygon");
      info(str(GEOSGeomTypeId(B)));
    }
    const GEOSGeometry *rA = GEOSGetExteriorRing(A);
    const GEOSGeometry *rB = GEOSGetExteriorRing(B);
    const GEOSCoordSequence *vA = GEOSGeom_getCoordSeq(rA);
    const GEOSCoordSequence *vB = GEOSGeom_getCoordSeq(rB);
    const int nA = GEOSGetNumCoordinates(rA);
    const int nB = GEOSGetNumCoordinates(rB);
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

        // Use either projection or closest vertex
        if (inside)
        {
          const double d = sqrt(Geometry::Dot2D(q - p, q - p));
          if (d < tol)
          {
            projections.push_back(p);
            projections.push_back(q);
          }
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
};

} // namespace DTCC_BUILDER

#endif
