// Copyright (C) 2022 Anders Logg
// Licensed under the MIT License

#ifndef GEOS_H
#define GEOS_H

#include <geos_c.h>

#include "Logging.h"
#include "Polygon.h"
#include "Timer.h"

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
class GEOS
{
public:
  /// Print GEOS geometry to terminal (in WKT format)
  static void Print(GEOSGeometry *geometry)
  {
    GEOSWKTWriter *writer = GEOSWKTWriter_create();
    GEOSWKTWriter_setTrim(writer, 1);
    char *out = GEOSWKTWriter_write(writer, geometry);
    Info(out);
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

    GEOSGeometry *ring = GEOSGeom_createLinearRing(sequence);
    GEOSGeometry *_polygon = GEOSGeom_createPolygon(ring, 0, 0);

    Print(ring);
    Print(_polygon);

    // Cleanup

    return 0;
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

    // FIXME: Figure out where to handle this
    initGEOS(geos_msg_handler, geos_msg_handler);

    GEOSGeometry *A = CreatePolygon(polygon0);

    // FIXME: Figure out where to handle this
    finishGEOS();

    Polygon p;
    return p;
  }
};

} // namespace DTCC

#endif
