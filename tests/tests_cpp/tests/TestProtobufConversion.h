#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

#include "dtcc.pb.h"
#include "protobuf/include/PointCloudMethods.h"

#include "PointCloud.h"
#include "Protobuf.h"

using namespace DTCC_BUILDER;

TEST_CASE("Protobuf to PointCloud")
{
  SECTION("Load protobuf pointcloud")
  {
    std::string pbFilePath = RootPath + "data/MinimalCase/pointcloud.las.pb";
    DTCC::PointCloud pb_pointCloud;
    std::ifstream input(pbFilePath, std::ios::in | std::ios::binary);
    std::string pbString(std::istreambuf_iterator<char>(input), {});

    REQUIRE(pb_pointCloud.ParseFromString(pbString));
    REQUIRE(pb_pointCloud.points().size() == 8148);
  }

  SECTION("Convert to PointCloud")
  {
    std::string pbFilePath = RootPath + "data/MinimalCase/pointcloud.las.pb";
    DTCC::PointCloud pb_pointCloud;
    std::fstream input(pbFilePath, std::ios::in | std::ios::binary);
    std::string pbString(std::istreambuf_iterator<char>(input), {});
    pb_pointCloud.ParseFromString(pbString);
    PointCloud pc = Protobuf::LoadPointCloud(pb_pointCloud);

    REQUIRE(pc.Points.size() == 8148);
    REQUIRE(pc.Classifications.size() == 8148);
  }
}
