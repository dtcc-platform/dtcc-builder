#include "Geometry.h"
#include "Point.h"
#include "Polyfix.h"
#include "Polygon.h"

using namespace DTCC;

TEST_CASE("Polygon Geometry")
{
  SECTION("Polygon Area")
  {
    Polygon polygon;
    polygon.Vertices.push_back(Point2D(0, 0));
    polygon.Vertices.push_back(Point2D(1, 0));
    polygon.Vertices.push_back(Point2D(1, 1));
    polygon.Vertices.push_back(Point2D(0, 1));
    REQUIRE(Geometry::PolygonArea(polygon) == Approx(1.0));

    polygon.Vertices.clear();
    polygon.Vertices.push_back(Point2D(0, 0));
    polygon.Vertices.push_back(Point2D(3, 0));
    polygon.Vertices.push_back(Point2D(3, 3));
    polygon.Vertices.push_back(Point2D(2, 3));
    polygon.Vertices.push_back(Point2D(2, 2));
    polygon.Vertices.push_back(Point2D(1, 2));
    polygon.Vertices.push_back(Point2D(1, 3));
    polygon.Vertices.push_back(Point2D(0, 3));
    REQUIRE(Geometry::PolygonArea(polygon) == Approx(8.0));
  }

  SECTION("Point in Polygon")
  {
    Polygon polygon;
    polygon.Vertices.push_back(Point2D(0, 0));
    polygon.Vertices.push_back(Point2D(3, 0));
    polygon.Vertices.push_back(Point2D(3, 3));
    polygon.Vertices.push_back(Point2D(2, 3));
    polygon.Vertices.push_back(Point2D(2, 2));
    polygon.Vertices.push_back(Point2D(1, 2));
    polygon.Vertices.push_back(Point2D(1, 3));
    polygon.Vertices.push_back(Point2D(0, 3));

    REQUIRE(Geometry::PolygonContains2D(polygon, Point2D(0, 0)));
    REQUIRE(Geometry::PolygonContains2D(polygon, Point2D(1, 1)));
    REQUIRE(!Geometry::PolygonContains2D(polygon, Point2D(1.5, 2.5)));
  }

  SECTION("Polygon Centroid")
  {
    Polygon polygon;
    polygon.Vertices.push_back(Point2D(0, 0));
    polygon.Vertices.push_back(Point2D(1, 0));
    polygon.Vertices.push_back(Point2D(1, 1));
    polygon.Vertices.push_back(Point2D(0, 1));
    auto centroid = Geometry::PolygonCenter2D(polygon);
    REQUIRE(centroid.x == 0.5);
    REQUIRE(centroid.y == 0.5);

    polygon.Vertices.clear();
    polygon.Vertices.push_back(Point2D(0, 0));
    polygon.Vertices.push_back(Point2D(3, 0));
    polygon.Vertices.push_back(Point2D(3, 3));
    polygon.Vertices.push_back(Point2D(2, 3));
    polygon.Vertices.push_back(Point2D(2, 2));
    polygon.Vertices.push_back(Point2D(1, 2));
    polygon.Vertices.push_back(Point2D(1, 3));
    polygon.Vertices.push_back(Point2D(0, 3));
    centroid = Geometry::PolygonCenter2D(polygon);
    REQUIRE(centroid.x == 1.5);
    REQUIRE(centroid.y == 2.0);
  }

  SECTION("Polygon Perimeter")
  {
    Polygon polygon;
    polygon.Vertices.push_back(Point2D(0, 0));
    polygon.Vertices.push_back(Point2D(1, 0));
    polygon.Vertices.push_back(Point2D(1, 1));
    polygon.Vertices.push_back(Point2D(0, 1));
    REQUIRE(Geometry::PolygonPerimeter2D(polygon) == Approx(4.0));

    polygon.Vertices.clear();
    polygon.Vertices.push_back(Point2D(0, 0));
    polygon.Vertices.push_back(Point2D(3, 0));
    polygon.Vertices.push_back(Point2D(3, 3));
    polygon.Vertices.push_back(Point2D(2, 3));
    polygon.Vertices.push_back(Point2D(2, 2));
    polygon.Vertices.push_back(Point2D(1, 2));
    polygon.Vertices.push_back(Point2D(1, 3));
    polygon.Vertices.push_back(Point2D(0, 3));

    REQUIRE(Geometry::PolygonPerimeter2D(polygon) == Approx(14.0));
  }

  SECTION("Subtract polygon origin")
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
}
