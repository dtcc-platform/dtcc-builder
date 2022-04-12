#include "GEOS.h"
#include "Polygon.h"

using namespace DTCC;

TEST_CASE("GEOS::MergePolygons")
{
  const double TOL = 0.1;
  GEOS::Init();

  // Note: These are the same tests as in the sandbox scripts
  // MergePolygons.py and MergePolygonsGEOS.py.

  SECTION("Test case 0")
  {
    Polygon p0;
    p0.Vertices.push_back(Point2D(0, 0));
    p0.Vertices.push_back(Point2D(1, 0));
    p0.Vertices.push_back(Point2D(1, 1));
    p0.Vertices.push_back(Point2D(0, 1));

    Polygon p1;
    for (const auto &p : p0.Vertices)
      p1.Vertices.push_back(p + Vector2D(1.1, 0.5));

    GEOS::MergePolygons(p0, p1, TOL);
  }

  SECTION("Test case 1")
  {
    Polygon p0;
    p0.Vertices.push_back(Point2D(0, 0));
    p0.Vertices.push_back(Point2D(1, 0));
    p0.Vertices.push_back(Point2D(1, 1));
    p0.Vertices.push_back(Point2D(0, 1));

    Polygon p1;
    for (const auto &p : p0.Vertices)
      p1.Vertices.push_back(p + Vector2D(0.6, 0.5));

    GEOS::MergePolygons(p0, p1, TOL);
  }
}
