#include "BoundingBox.h"
#include "PointCloudProcessor.h"
#include "model/Point.h"
#include "model/PointCloud.h"

using namespace DTCC_BUILDER;

TEST_CASE("POINT_CLOUD")
{

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

  SECTION("Used Classifications")
  {
    PointCloud pc;
    pc.Points.push_back(Vector3D(0, 0, 0));
    pc.Classifications.push_back(0);
    pc.Points.push_back(Vector3D(1, 0, 0));
    pc.Classifications.push_back(1);
    pc.Points.push_back(Vector3D(1, 0, 0));
    pc.Classifications.push_back(1);
    pc.Points.push_back(Vector3D(1, 0, 0));
    pc.Classifications.push_back(1);
    pc.Points.push_back(Vector3D(2, 0, 0));
    pc.Classifications.push_back(2);
    pc.Points.push_back(Vector3D(2, 0, 0));
    pc.Classifications.push_back(2);
    pc.Points.push_back(Vector3D(2, 0, 0));
    pc.Classifications.push_back(2);

    pc.BuildHasClassifications();

    REQUIRE(pc.HasClassification(0));
    REQUIRE(pc.HasClassification(1));
    REQUIRE(pc.HasClassification(2));
    REQUIRE(!pc.HasClassification(3));
  }
}

TEST_CASE("Outlier remover")
{
  PointCloud pc;
  pc.Points.push_back(Vector3D(0, 0, 0));
  pc.Points.push_back(Vector3D(0.5, 0.5, 0));
  pc.Points.push_back(Vector3D(0.5, 0.5, 0.5));
  pc.Points.push_back(Vector3D(1, 1, 1));
  pc.Points.push_back(Vector3D(1.5, 1.5, 1));
  pc.Points.push_back(Vector3D(1.5, 1.5, 1.5));
  pc.Points.push_back(Vector3D(10, 10, 10));
  SECTION("Nearest Neighbours")
  {
    auto knn = PointCloudProcessor::KNNNearestNeighbours(pc.Points, 3);
    REQUIRE(knn.at(0).size() == 3);
    REQUIRE(knn.at(5).size() == 3);
    REQUIRE(knn.at(0).at(0) == std::sqrt(0.5 * 0.5 + 0.5 * 0.5));
    REQUIRE(knn.at(0).at(1) == std::sqrt(0.5 * 0.5 + 0.5 * 0.5 + 0.5 * 0.5));
    REQUIRE(knn.at(6).at(0) ==
            std::sqrt((10 - 1.5) * (10 - 1.5) + (10 - 1.5) * (10 - 1.5) +
                      (10 - 1.5) * (10 - 1.5)));
    REQUIRE(knn.at(6).at(1) ==
            std::sqrt((10 - 1.5) * (10 - 1.5) + (10 - 1.5) * (10 - 1.5) +
                      (10 - 1) * (10 - 1)));
  }

  SECTION("Outlier Remover")
  {
    auto outliers =
        PointCloudProcessor::StatisticalOutlierFinder(pc.Points, 3, 1.5);
    REQUIRE(outliers.size() == 1);
    REQUIRE(outliers[0] == 6);

    REQUIRE(pc.Points.size() == 7);
    PointCloudProcessor::StatisticalOutlierRemover(pc, 3, 1.5);
    REQUIRE(pc.Points.size() == 6);
  }

  SECTION("Parse Scan Flag")
  {
    auto flag1 = PointCloudProcessor::parseScanFlag(9);
    auto flag2 = PointCloudProcessor::parseScanFlag(26);
    REQUIRE(flag1.first == 1);
    REQUIRE(flag1.second == 1);
    REQUIRE(flag2.first == 2);
    REQUIRE(flag2.second == 3);
  }

  SECTION("Pack Scan Flag")
  {
    auto flag1 = PointCloudProcessor::packScanFlag(1, 1);
    auto flag2 = PointCloudProcessor::packScanFlag(2, 3);
    REQUIRE(flag1 == 9);

    auto parse1 = PointCloudProcessor::parseScanFlag(flag1);
    auto parse2 = PointCloudProcessor::parseScanFlag(flag2);
    REQUIRE(parse1.first == 1);
    REQUIRE(parse1.second == 1);
    REQUIRE(parse2.first == 2);
    REQUIRE(parse2.second == 3);
  }
}

TEST_CASE("Vegetation filter")
{
  PointCloud pc;
  pc.Points.push_back(Vector3D(0, 0, 0));
  pc.Classifications.push_back(1);
  pc.Points.push_back(Vector3D(0, 0, 1));
  pc.Classifications.push_back(1);
  pc.Points.push_back(Vector3D(0, 0, 2));
  pc.Classifications.push_back(1);
  pc.Points.push_back(Vector3D(0, 0, 3));
  pc.Classifications.push_back(1);

  SECTION("No flags")
  {
    // No flags, do nothing
    size_t pre_filter = pc.Points.size();
    pc = PointCloudProcessor::remove_vegetation(pc);
    REQUIRE(pc.Points.size() == pre_filter);
  }

  SECTION("Filter")
  {
    pc.ScanFlags.push_back(PointCloudProcessor::packScanFlag(1, 1));
    pc.ScanFlags.push_back(PointCloudProcessor::packScanFlag(1, 2));
    pc.ScanFlags.push_back(PointCloudProcessor::packScanFlag(2, 2));
    pc.ScanFlags.push_back(PointCloudProcessor::packScanFlag(2, 3));

    pc = PointCloudProcessor::remove_vegetation(pc);
    REQUIRE(pc.Points.size() == 2);
  }
}

TEST_CASE("RANSAC filter")
{
  std::vector<Point3D> points;
  points.push_back(Vector3D(0, 0, 25));
  points.push_back(Vector3D(0, 0, -25));
  for (int x = 0; x < 10; x++)
  {
    for (int y = 0; y < 10; y++)
    {
      points.push_back(Vector3D(x, y, 0));
    }
  }
  SECTION("Outliers")
  {
    auto outliers = PointCloudProcessor::RANSAC_OutlierFinder(points, 0.1, 150);
    REQUIRE(outliers.size() == 2);
    REQUIRE(outliers[0] == 0);
    REQUIRE(outliers[1] == 1);
  }
}
