#include "Polygon.h"
#include "SHP.h"
#include "Utils.h"

using namespace DTCC_BUILDER;

TEST_CASE("Create own uuids")
{
  // Create regex matching UUIDs
  std::string hx = "[[:xdigit:]]";
  std::string uuid_regex(hx + "{8}-" + hx + "{4}-" + hx + "{4}-" + hx + "{4}-" +
                         hx + "{12}");

  // Test utility function
  std::string uuid = Utils::CreateUUID();
  REQUIRE(std::regex_match(uuid, std::regex(uuid_regex)));

  // Test property map sample files
  std::vector<Polygon> polygons;
  std::vector<std::string> uuids;
  std::vector<int> entity_ids;
  std::string filename =
      root_path + "data/create-uuid-propertymap/PropertyMap.shp";
  SHP::Read(polygons, filename, &uuids, &entity_ids);

  REQUIRE(uuids[0] == "0d1a05a7-c188-42ab-adad-35b359f1b920");
  REQUIRE(
      std::regex_match(uuids[1], std::regex(uuid_regex + "::DTCCGenerated")));
}

TEST_CASE("Add UUID-N on duplicate uuids")
{
  std::string filename = root_path + "data/duplicate-shp-uuids/by_14.shp";
  std::vector<Polygon> polygons;
  std::vector<std::string> uuids;
  std::vector<int> entityID;

  SECTION("Handle duplicate uuids")
  {
    SHP::Read(polygons, filename, &uuids, &entityID);
    REQUIRE(uuids[118] == "e5e21821-9108-4762-9caa-8e612f81febb-1");
    REQUIRE(uuids[150] == "e5e21821-9108-4762-9caa-8e612f81febb-2");
    // Only unique UUIDs
    REQUIRE(uuids.size() ==
            std::set<std::string>(uuids.begin(), uuids.end()).size());
  }

  SECTION("Allow duplicate uuids")
  {
    const int def_uuid_length = 36;
    SHP::Read(polygons, filename, &uuids, &entityID, nullptr, false);
    // All UUIDs of default length
    REQUIRE(std::count_if(uuids.begin(), uuids.end(),
                          [](const std::string &UUID) {
                            return UUID.length() == def_uuid_length;
                          }) == uuids.size());
  }
}