#include "Point.h"
#include "Polyfix.h"
#include "Polygon.h"

using namespace DTCC;

TEST_CASE("Subtract polygon origin")
{
  std::vector<Point2D> vertices = {Point2D(3, 7), Point2D(2, 4.5)};
  Polygon p;
  p.Vertices = vertices;
  Polyfix::Transform(p, Point2D(1, 2));
  Polyfix::Transform(vertices, Point2D(1, 2));
  for (const auto &vertices2 : {vertices, p.Vertices})
  {
    REQUIRE(vertices2[0].x == 2);
    REQUIRE(vertices2[0].y == 5);
    REQUIRE(vertices2[1].x == 1);
    REQUIRE(vertices2[1].y == 2.5);
  }
}