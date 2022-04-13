#include "GEOS.h"
#include "Polygon.h"

using namespace DTCC;

// Check sum for comparing against sandbox/MergePolygonsGEOS.py
double CheckSum(const Polygon &polygon)
{
  double s = 0.0;
  for (const auto &p : polygon.Vertices)
    s += p.x + p.y;
  return s;
}

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

    Polygon pm = GEOS::MergePolygons(p0, p1, TOL);
    REQUIRE(CheckSum(pm) == Approx(14.3));
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

    Polygon pm = GEOS::MergePolygons(p0, p1, TOL);
    REQUIRE(CheckSum(pm) == Approx(12.4));
  }
}

/*
Checksum: 14.3
Checksum: 12.4
Checksum: 14.109756
Checksum: 6.0
Checksum: 4873.52996
Checksum: 14.344555
Checksum: 17905.235323
Checksum: 4489.04597
Checksum: 2326.688678
Checksum: 9424.87874
Checksum: 13163.017939
Checksum: 16962.607815000003
Checksum: 7542.28551
Checksum: 16767.20206
*/
