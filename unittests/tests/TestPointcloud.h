#include "BoundingBox.h"
#include "CSV.h"
#include "LAS.h"
#include "Point.h"
#include "PointCloud.h"
#include "PointCloudProcessor.h"

using namespace DTCC;

TEST_CASE("POINT_CLOUD")
{
  SECTION("READ LAS")
  {
    PointCloud pc;
    LAS::Read(pc, RootPath + "data/minimal_las.las");

    REQUIRE(pc.Points.size() == 10);
    for (size_t i = 0; i < pc.Points.size(); i++)
    {
      REQUIRE(pc.Classifications[i] == Approx(pc.Points[i].x).margin(1e-6));
    }
  }

  SECTION("BOUNDS")
  {
    BoundingBox2D bb;
    LAS::Bounds(bb, RootPath + "data/minimal_las.las");
    REQUIRE(bb.P.x == 0);
    REQUIRE(bb.Q.x == 9);
  }

  SECTION("ClassificationFilter")
  {
    PointCloud pc;
    pc.Points.push_back(Vector3D(0, 0, 0));
    pc.Classifications.push_back(0);
    pc.Points.push_back(Vector3D(1, 0, 0));
    pc.Classifications.push_back(1);
    pc.Points.push_back(Vector3D(2, 0, 0));
    pc.Classifications.push_back(2);

    PointCloud out_pc = PointCloudProcessor::ClassificationFilter(pc, {1, 2});
    REQUIRE(out_pc.Points.size() == 2);
    REQUIRE(out_pc.Classifications.size() == 2);
    REQUIRE(out_pc.Points[0].x == 1);
    REQUIRE(out_pc.Points[1].x == 2);
  }
}

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

TEST_CASE("Read from CSV instead of LAS/LAZ")
{
  std::string filename =
      RootPath + "data/read-from-csv-instead-of-laz/PointCloudTest.csv";
  PointCloud pointCloud;
  CSV::Read(pointCloud, filename);

  SECTION("PointCloud vertices")
  {
    Point3D v1 = pointCloud.Points[0];
    REQUIRE(v1.x == 317228.73);
    REQUIRE(v1.y == 6397500.00);
    REQUIRE(v1.z == 26.16);
  }

  SECTION("PointCloud colors")
  {
    Color c1 = pointCloud.Colors[0];
    for (double ch : {c1.R, c1.G, c1.B})
      REQUIRE(ch == 0);
  }

  SECTION("PointCloud classification")
  {
    auto classification = pointCloud.Classifications[0];
    REQUIRE(classification == 1);
  }

  SECTION("Read points only within bounding box")
  {
    pointCloud.Clear();
    BoundingBox2D bbox(Point2D(315500, 6397510), Point2D(317000, 6399000));
    CSV::Read(pointCloud, filename, bbox);
    for (const auto &p : pointCloud.Points)
    {
      REQUIRE(p.x >= bbox.P.x);
      REQUIRE(p.y >= bbox.P.y);
      REQUIRE(p.x <= bbox.Q.x);
      REQUIRE(p.y <= bbox.Q.y);
    }
  }

  SECTION("Read only points of certain classification")
  {
    pointCloud.Clear();
    std::vector<int> groundWaterPts{2, 9};
    CSV::Read(pointCloud, filename, groundWaterPts);
    // 1 is only other present classification
    for (const auto &c : pointCloud.Classifications)
      REQUIRE(c != 1);
  }
}