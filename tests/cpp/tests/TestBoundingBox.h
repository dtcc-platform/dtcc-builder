#include "BoundingBox.h"

using namespace DTCC_BUILDER;

TEST_CASE("BoundingBox2D")
{
  Polygon p;
  p.vertices = {Vector2D(1, 2), Vector2D(0, 3), Vector2D(2, 1)};
  BoundingBox2D bbox_p(p);
  BoundingBox2D bbox_v(p.vertices);

  for (const auto &bbox : {bbox_p, bbox_v})
  {
    REQUIRE(bbox.P.x == 0);
    REQUIRE(bbox.P.y == 1);
    REQUIRE(bbox.Q.x == 2);
    REQUIRE(bbox.Q.y == 3);
  }

  Vector2D p1 = Vector2D(0, 0);
  Vector2D p2 = Vector2D(5, 5);
  Vector2D p3 = Vector2D(-5, -5);
  Vector2D p4 = Vector2D(-1, -1);
  Vector2D p5 = Vector2D(-2, -2);
  Vector2D p6 = Vector2D(2, 2);
  Vector2D p7 = Vector2D(1, 1);

  SECTION("UNION")
  {
    BoundingBox2D bb1 = BoundingBox2D(p1, p2);
    BoundingBox2D bb2 = BoundingBox2D(p3, p4);
    bb1.union_with(bb2);
    REQUIRE(bb1.P.x == -5);
    REQUIRE(bb1.P.y == -5);
    REQUIRE(bb1.Q.x == 5);
    REQUIRE(bb1.Q.y == 5);
  }

  SECTION("INTERSECT")
  {
    BoundingBox2D bb1 = BoundingBox2D(p1, p2);
    BoundingBox2D bb2 = BoundingBox2D(p3, p4);
    bb1.intersect(bb2);
    REQUIRE(bb1.area() == 0);

    BoundingBox2D bb3 = BoundingBox2D(p1, p2);
    BoundingBox2D bb4 = BoundingBox2D(p5, p6);
    bb3.intersect(bb4);
    REQUIRE(bb3.P.x == 0);
    REQUIRE(bb3.P.y == 0);
    REQUIRE(bb3.Q.x == 2);
    REQUIRE(bb3.Q.y == 2);

    BoundingBox2D bb5 = BoundingBox2D(p1, p2);
    BoundingBox2D bb6 = BoundingBox2D(p7, p6);
    bb5.intersect(bb6);
    REQUIRE(bb5.P.x == bb6.P.x);
    REQUIRE(bb5.P.y == bb6.P.y);
    REQUIRE(bb5.Q.x == bb6.Q.x);
    REQUIRE(bb5.Q.y == bb6.Q.y);
  }
}