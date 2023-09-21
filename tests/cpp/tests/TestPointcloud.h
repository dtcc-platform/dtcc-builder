#include "BoundingBox.h"
#include "PointCloudProcessor.h"
#include "model/PointCloud.h"
#include "model/Vector.h"

using namespace DTCC_BUILDER;

TEST_CASE("POINT_CLOUD")
{

  SECTION("classification_filter")
  {
    PointCloud pc;
    pc.points.push_back(Vector3D(0, 0, 0));
    pc.classifications.push_back(0);
    pc.points.push_back(Vector3D(1, 0, 0));
    pc.classifications.push_back(1);
    pc.points.push_back(Vector3D(2, 0, 0));
    pc.classifications.push_back(2);

    PointCloud out_pc = PointCloudProcessor::classification_filter(pc, {1, 2});
    REQUIRE(out_pc.points.size() == 2);
    REQUIRE(out_pc.classifications.size() == 2);
    REQUIRE(out_pc.points[0].x == 1);
    REQUIRE(out_pc.points[1].x == 2);
  }

  SECTION("Used classifications")
  {
    PointCloud pc;
    pc.points.push_back(Vector3D(0, 0, 0));
    pc.classifications.push_back(0);
    pc.points.push_back(Vector3D(1, 0, 0));
    pc.classifications.push_back(1);
    pc.points.push_back(Vector3D(1, 0, 0));
    pc.classifications.push_back(1);
    pc.points.push_back(Vector3D(1, 0, 0));
    pc.classifications.push_back(1);
    pc.points.push_back(Vector3D(2, 0, 0));
    pc.classifications.push_back(2);
    pc.points.push_back(Vector3D(2, 0, 0));
    pc.classifications.push_back(2);
    pc.points.push_back(Vector3D(2, 0, 0));
    pc.classifications.push_back(2);

    pc.build_has_classifications();

    REQUIRE(pc.has_classification(0));
    REQUIRE(pc.has_classification(1));
    REQUIRE(pc.has_classification(2));
    REQUIRE(!pc.has_classification(3));
  }
}

TEST_CASE("Outlier remover")
{
  PointCloud pc;
  pc.points.push_back(Vector3D(0, 0, 0));
  pc.points.push_back(Vector3D(0.5, 0.5, 0));
  pc.points.push_back(Vector3D(0.5, 0.5, 0.5));
  pc.points.push_back(Vector3D(1, 1, 1));
  pc.points.push_back(Vector3D(1.5, 1.5, 1));
  pc.points.push_back(Vector3D(1.5, 1.5, 1.5));
  pc.points.push_back(Vector3D(10, 10, 10));
  SECTION("nearest Neighbours")
  {
    auto knn = PointCloudProcessor::knn_nearest_neighbours(pc.points, 3);
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
        PointCloudProcessor::statistical_outlier_finder(pc.points, 3, 1.5);
    REQUIRE(outliers.size() == 1);
    REQUIRE(outliers[0] == 6);

    REQUIRE(pc.points.size() == 7);
    PointCloudProcessor::statistical_outlier_remover(pc, 3, 1.5);
    REQUIRE(pc.points.size() == 6);
  }

  SECTION("Parse Scan Flag")
  {
    auto flag1 = PointCloudProcessor::parse_scan_flag(9);
    auto flag2 = PointCloudProcessor::parse_scan_flag(26);
    REQUIRE(flag1.first == 1);
    REQUIRE(flag1.second == 1);
    REQUIRE(flag2.first == 2);
    REQUIRE(flag2.second == 3);
  }

  SECTION("Pack Scan Flag")
  {
    auto flag1 = PointCloudProcessor::pack_scan_flag(1, 1);
    auto flag2 = PointCloudProcessor::pack_scan_flag(2, 3);
    REQUIRE(flag1 == 9);

    auto parse1 = PointCloudProcessor::parse_scan_flag(flag1);
    auto parse2 = PointCloudProcessor::parse_scan_flag(flag2);
    REQUIRE(parse1.first == 1);
    REQUIRE(parse1.second == 1);
    REQUIRE(parse2.first == 2);
    REQUIRE(parse2.second == 3);
  }
}

TEST_CASE("Vegetation filter")
{
  PointCloud pc;
  pc.points.push_back(Vector3D(0, 0, 0));
  pc.classifications.push_back(1);
  pc.points.push_back(Vector3D(0, 0, 1));
  pc.classifications.push_back(1);
  pc.points.push_back(Vector3D(0, 0, 2));
  pc.classifications.push_back(1);
  pc.points.push_back(Vector3D(0, 0, 3));
  pc.classifications.push_back(1);

  SECTION("No flags")
  {
    // No flags, do nothing
    size_t pre_filter = pc.points.size();
    pc = PointCloudProcessor::remove_vegetation(pc);
    REQUIRE(pc.points.size() == pre_filter);
  }

  SECTION("Filter")
  {
    pc.scan_flags.push_back(PointCloudProcessor::pack_scan_flag(1, 1));
    pc.scan_flags.push_back(PointCloudProcessor::pack_scan_flag(1, 2));
    pc.scan_flags.push_back(PointCloudProcessor::pack_scan_flag(2, 2));
    pc.scan_flags.push_back(PointCloudProcessor::pack_scan_flag(2, 3));

    pc = PointCloudProcessor::remove_vegetation(pc);
    REQUIRE(pc.points.size() == 2);
  }
}

TEST_CASE("RANSAC filter")
{
  std::vector<Vector3D> points;
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
    auto outliers =
        PointCloudProcessor::ransac_outlier_finder(points, 0.1, 150);
    REQUIRE(outliers.size() == 2);
    REQUIRE(outliers[0] == 0);
    REQUIRE(outliers[1] == 1);
  }
}
