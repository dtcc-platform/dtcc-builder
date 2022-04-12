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

  /// Create GEOS polygon from Polygon
  static GEOSGeometry *CreatePolygon(const Polygon &polygon)
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
    GEOSGeometry *_polygon = GEOSGeom_createPolygon(ring, 0, 0);

    return _polygon;
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

    // Create GEOS polygons
    GEOSGeometry *A = CreatePolygon(polygon0);
    GEOSGeometry *B = CreatePolygon(polygon1);

    // Simplify geometries
    GEOSGeometry *_A = GEOSSimplify(A, tol);
    GEOSGeometry *_B = GEOSSimplify(B, tol);

    // Compute union
    GEOSGeometry *C = GEOSUnionPrec(A, B, EPS);
    GEOSGeometry *_C = GEOSSimplify(C, tol);

    // Cleanup
    GEOSGeom_destroy(A);
    GEOSGeom_destroy(B);
    GEOSGeom_destroy(C);
    GEOSGeom_destroy(_A);
    GEOSGeom_destroy(_B);
    GEOSGeom_destroy(_C);

    Polygon p;
    return p;
  }
  };

  } // namespace DTCC

#endif
