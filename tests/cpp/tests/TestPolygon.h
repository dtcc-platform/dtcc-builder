#include "Geometry.h"
#include "Polyfix.h"
#include "model/Polygon.h"
#include "model/Vector.h"

using namespace DTCC_BUILDER;

TEST_CASE("Polygon Geometry")
{
  SECTION("Polygon area")
  {
    Polygon polygon;
    polygon.vertices.push_back(Vector2D(0, 0));
    polygon.vertices.push_back(Vector2D(1, 0));
    polygon.vertices.push_back(Vector2D(1, 1));
    polygon.vertices.push_back(Vector2D(0, 1));
    REQUIRE(Geometry::polygon_area(polygon) == Approx(1.0));

    polygon.vertices.clear();
    polygon.vertices.push_back(Vector2D(0, 0));
    polygon.vertices.push_back(Vector2D(3, 0));
    polygon.vertices.push_back(Vector2D(3, 3));
    polygon.vertices.push_back(Vector2D(2, 3));
    polygon.vertices.push_back(Vector2D(2, 2));
    polygon.vertices.push_back(Vector2D(1, 2));
    polygon.vertices.push_back(Vector2D(1, 3));
    polygon.vertices.push_back(Vector2D(0, 3));
    REQUIRE(Geometry::polygon_area(polygon) == Approx(8.0));
  }

  SECTION("Point in Polygon")
  {
    Polygon polygon;
    polygon.vertices.push_back(Vector2D(0, 0));
    polygon.vertices.push_back(Vector2D(3, 0));
    polygon.vertices.push_back(Vector2D(3, 3));
    polygon.vertices.push_back(Vector2D(2, 3));
    polygon.vertices.push_back(Vector2D(2, 2));
    polygon.vertices.push_back(Vector2D(1, 2));
    polygon.vertices.push_back(Vector2D(1, 3));
    polygon.vertices.push_back(Vector2D(0, 3));

    REQUIRE(Geometry::polygon_contains_2d(polygon, Vector2D(0, 0)));
    REQUIRE(Geometry::polygon_contains_2d(polygon, Vector2D(1, 1)));
    REQUIRE(!Geometry::polygon_contains_2d(polygon, Vector2D(1.5, 2.5)));
  }

  SECTION("Polygon Centroid")
  {
    Polygon polygon;
    polygon.vertices.push_back(Vector2D(0, 0));
    polygon.vertices.push_back(Vector2D(1, 0));
    polygon.vertices.push_back(Vector2D(1, 1));
    polygon.vertices.push_back(Vector2D(0, 1));
    auto centroid = Geometry::polygon_center_2d(polygon);
    REQUIRE(centroid.x == 0.5);
    REQUIRE(centroid.y == 0.5);

    polygon.vertices.clear();
    polygon.vertices.push_back(Vector2D(0, 0));
    polygon.vertices.push_back(Vector2D(3, 0));
    polygon.vertices.push_back(Vector2D(3, 3));
    polygon.vertices.push_back(Vector2D(2, 3));
    polygon.vertices.push_back(Vector2D(2, 2));
    polygon.vertices.push_back(Vector2D(1, 2));
    polygon.vertices.push_back(Vector2D(1, 3));
    polygon.vertices.push_back(Vector2D(0, 3));
    centroid = Geometry::polygon_center_2d(polygon);
    REQUIRE(centroid.x == 1.5);
    REQUIRE(centroid.y == 2.0);
  }

  SECTION("Polygon Perimeter")
  {
    Polygon polygon;
    polygon.vertices.push_back(Vector2D(0, 0));
    polygon.vertices.push_back(Vector2D(1, 0));
    polygon.vertices.push_back(Vector2D(1, 1));
    polygon.vertices.push_back(Vector2D(0, 1));
    REQUIRE(Geometry::polygon_perimeter_2d(polygon) == Approx(4.0));

    polygon.vertices.clear();
    polygon.vertices.push_back(Vector2D(0, 0));
    polygon.vertices.push_back(Vector2D(3, 0));
    polygon.vertices.push_back(Vector2D(3, 3));
    polygon.vertices.push_back(Vector2D(2, 3));
    polygon.vertices.push_back(Vector2D(2, 2));
    polygon.vertices.push_back(Vector2D(1, 2));
    polygon.vertices.push_back(Vector2D(1, 3));
    polygon.vertices.push_back(Vector2D(0, 3));

    REQUIRE(Geometry::polygon_perimeter_2d(polygon) == Approx(14.0));
  }

  SECTION("Subtract polygon origin")
  {
    std::vector<Vector2D> vertices = {Vector2D(3, 7), Vector2D(2, 4.5)};
    Polygon p;
    p.vertices = vertices;
    Polyfix::transform(p, Vector2D(1, 2));
    Polyfix::transform(vertices, Vector2D(1, 2));
    for (const auto &vertices2 : {vertices, p.vertices})
    {
      REQUIRE(vertices2[0].x == 2);
      REQUIRE(vertices2[0].y == 5);
      REQUIRE(vertices2[1].x == 1);
      REQUIRE(vertices2[1].y == 2.5);
    }
  }

  SECTION("Polygon Intersects")
  {
    Polygon polygon;
    polygon.vertices.push_back(Vector2D(0, 0));
    polygon.vertices.push_back(Vector2D(1, 0));
    polygon.vertices.push_back(Vector2D(1, 1));
    polygon.vertices.push_back(Vector2D(0, 1));

    Polygon polygon2;
    polygon2.vertices.push_back(Vector2D(0.5, 0.5));
    polygon2.vertices.push_back(Vector2D(1.5, 0.5));
    polygon2.vertices.push_back(Vector2D(1.5, 1.5));
    polygon2.vertices.push_back(Vector2D(0.5, 1.5));

    Polygon polygon3;
    polygon3.vertices.push_back(Vector2D(-10, -10));
    polygon3.vertices.push_back(Vector2D(10, -10));
    polygon3.vertices.push_back(Vector2D(10, 10));
    polygon3.vertices.push_back(Vector2D(-10, 10));

    Polygon polygon4;
    polygon4.vertices.push_back(Vector2D(2, 2));
    polygon4.vertices.push_back(Vector2D(3, 2));
    polygon4.vertices.push_back(Vector2D(3, 3));
    polygon4.vertices.push_back(Vector2D(2, 3));

    REQUIRE(Geometry::intersects_2d(polygon, polygon));
    REQUIRE(Geometry::intersects_2d(polygon, polygon2));
    REQUIRE(Geometry::intersects_2d(polygon, polygon3));
    REQUIRE(!Geometry::intersects_2d(polygon, polygon4));
  }
}
