#include "Point.h"
#include "RoofSegmentation.h"

using namespace DTCC_BUILDER;

TEST_CASE("RoofSegmentation")
{
  std::vector<Point3D> points;
  for (int x = 0; x < 10; x++)
  {
    for (int y = 0; y < 10; y++)
    {
      points.push_back(Point3D(x, y, y * 0.5));
    }
  }
  for (int x = 0; x < 10; x++)
  {
    for (int y = 0; y < 10; y++)
    {
      points.push_back(Point3D(x + 10, y, 10 - y * 0.5));
    }
  }
  SECTION("Handle no data")
  {
    std::vector<Point3D> points;
    auto regions = RoofSegmentation::RegionGrowingSegmentation(points, 4);
    REQUIRE(regions.size() == 0);
  }
  SECTION("Region growing segmentation")
  {
    auto regions = RoofSegmentation::RegionGrowingSegmentation(points, 4);
    REQUIRE(regions.size() == 2);
    REQUIRE(regions[0].size() > 95);
    REQUIRE(regions[1].size() > 95);
    for (auto idx : regions[0])
    {
      REQUIRE(idx < 100);
    }
    for (auto idx : regions[1])
    {
      REQUIRE(idx >= 100);
    }
  }

  SECTION("Optimal Angle")
  {
    double angle =
        RoofSegmentation::SearchOptimateRegionSegmentationAngle(points, 4);
    REQUIRE(angle >= 0.5);
    REQUIRE(angle < 1.0);
  }
}