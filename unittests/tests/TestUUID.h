#include "Polygon.h"
#include "SHP.h"
#include "Utils.h"

using namespace DTCC_BUILDER;

TEST_CASE("Create own UUIDs")
{
  // Create regex matching UUIDs
  std::string hx = "[[:xdigit:]]";
  std::string uuidRegex(hx + "{8}-" + hx + "{4}-" + hx + "{4}-" + hx + "{4}-" +
                        hx + "{12}");

  // Test utility function
  std::string uuid = Utils::CreateUUID();
  REQUIRE(std::regex_match(uuid, std::regex(uuidRegex)));

  // Test property map sample files
  std::vector<Polygon> polygons;
  std::vector<std::string> uuids;
  std::vector<int> entityIDs;
  std::string filename =
      RootPath + "data/create-uuid-propertymap/PropertyMap.shp";
  SHP::Read(polygons, filename, &uuids, &entityIDs);

  REQUIRE(uuids[0] == "0d1a05a7-c188-42ab-adad-35b359f1b920");
  REQUIRE(
      std::regex_match(uuids[1], std::regex(uuidRegex + "::DTCCGenerated")));
}

TEST_CASE("Add UUID-N on duplicate UUIDs")
{
  std::string fileName = RootPath + "data/duplicate-shp-uuids/by_14.shp";
  std::vector<Polygon> polygons;
  std::vector<std::string> UUIDs;
  std::vector<int> entityID;

  SECTION("Handle duplicate UUIDs")
  {
    SHP::Read(polygons, fileName, &UUIDs, &entityID);
    REQUIRE(UUIDs[118] == "e5e21821-9108-4762-9caa-8e612f81febb-1");
    REQUIRE(UUIDs[150] == "e5e21821-9108-4762-9caa-8e612f81febb-2");
    // Only unique UUIDs
    REQUIRE(UUIDs.size() ==
            std::set<std::string>(UUIDs.begin(), UUIDs.end()).size());
  }

  SECTION("Allow duplicate UUIDs")
  {
    const int defUUIDLength = 36;
    SHP::Read(polygons, fileName, &UUIDs, &entityID, nullptr, false);
    // All UUIDs of default length
    REQUIRE(std::count_if(UUIDs.begin(), UUIDs.end(),
                          [](const std::string &UUID) {
                            return UUID.length() == defUUIDLength;
                          }) == UUIDs.size());
  }
}