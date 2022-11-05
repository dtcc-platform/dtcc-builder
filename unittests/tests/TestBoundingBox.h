#include "BoundingBox.h"

using namespace DTCC_BUILDER;

TEST_CASE("BoundingBox2D")
{
  Polygon p;
  p.Vertices = {Point2D(1, 2), Point2D(0, 3), Point2D(2, 1)};
  BoundingBox2D bboxP(p);
  BoundingBox2D bboxV(p.Vertices);

  for (const auto &bbox : {bboxP, bboxV})
  {
    REQUIRE(bbox.P.x == 0);
    REQUIRE(bbox.P.y == 1);
    REQUIRE(bbox.Q.x == 2);
    REQUIRE(bbox.Q.y == 3);
  }

  Point2D p1 = Point2D(0, 0);
  Point2D p2 = Point2D(5, 5);
  Point2D p3 = Point2D(-5, -5);
  Point2D p4 = Point2D(-1, -1);
  Point2D p5 = Point2D(-2, -2);
  Point2D p6 = Point2D(2, 2);
  Point2D p7 = Point2D(1, 1);

  SECTION("UNION")
  {
    BoundingBox2D bb1 = BoundingBox2D(p1, p2);
    BoundingBox2D bb2 = BoundingBox2D(p3, p4);
    bb1.Union(bb2);
    REQUIRE(bb1.P.x == -5);
    REQUIRE(bb1.P.y == -5);
    REQUIRE(bb1.Q.x == 5);
    REQUIRE(bb1.Q.y == 5);
  }

  SECTION("INTERSECT")
  {
    BoundingBox2D bb1 = BoundingBox2D(p1, p2);
    BoundingBox2D bb2 = BoundingBox2D(p3, p4);
    bb1.Intersect(bb2);
    REQUIRE(bb1.Area() == 0);

    BoundingBox2D bb3 = BoundingBox2D(p1, p2);
    BoundingBox2D bb4 = BoundingBox2D(p5, p6);
    bb3.Intersect(bb4);
    REQUIRE(bb3.P.x == 0);
    REQUIRE(bb3.P.y == 0);
    REQUIRE(bb3.Q.x == 2);
    REQUIRE(bb3.Q.y == 2);

    BoundingBox2D bb5 = BoundingBox2D(p1, p2);
    BoundingBox2D bb6 = BoundingBox2D(p7, p6);
    bb5.Intersect(bb6);
    REQUIRE(bb5.P.x == bb6.P.x);
    REQUIRE(bb5.P.y == bb6.P.y);
    REQUIRE(bb5.Q.x == bb6.Q.x);
    REQUIRE(bb5.Q.y == bb6.Q.y);
  }
}