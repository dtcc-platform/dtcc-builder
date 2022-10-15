#include "Polygon.h"
#include "SHP.h"

using namespace DTCCBUILDER;

TEST_CASE("Read SHP")
{
  SECTION("Load polygons")
  {
    std::vector<Polygon> buildings;
    SHP::Read(buildings, RootPath + "data/buildings/test_buildings.shp");
    REQUIRE(buildings.size() == 5);
  }
}