#include "BaseArea.h"
#include "District.h"
#include "JSON.h"

using namespace DTCC_BUILDER;

TEST_CASE("datamodel")
{
  const std::string fileName1 = root_path + "data/CityExample.json";
  const std::string fileName2 = root_path + "data/CityExample2.json";

  for (bool firstRun : {true, false})
  {
    District district;
    BaseArea baseArea;
    nlohmann::json json;

    JSON::Read(json, firstRun ? fileName1 : fileName2);

    if (firstRun)
    {
      JSON::Deserialize(district, json, "606");

      REQUIRE(district.AreaID == "606");
      REQUIRE(district.Name == "Hammarkullen");
      REQUIRE(district.Footprint.Vertices[0].x == Approx(-446.4952344278572));
      REQUIRE(district.Footprint.Vertices[0].y == Approx(150.96198354940861));
      REQUIRE(district.PrimaryAreas.size() == 1);

      PrimaryArea primaryArea = district.PrimaryAreas[0];
      REQUIRE(primaryArea.AreaID == "606");
      REQUIRE(primaryArea.Name == "Hammarkullen");
      REQUIRE(primaryArea.DistrictAreaID == "606");
      REQUIRE(primaryArea.Footprint.Vertices[0].x ==
              Approx(-446.4952344278572));
      REQUIRE(primaryArea.Footprint.Vertices[0].y ==
              Approx(150.96198354940861));

      baseArea = primaryArea.BaseAreas[0];
    }
    else
      JSON::Deserialize(baseArea, json, "60605");

    REQUIRE(baseArea.AreaID == "60605");
    REQUIRE(baseArea.PrimaryAreaID == "606");
    REQUIRE(baseArea.Footprint.Vertices[0].x == Approx(214.59035302355187));
    REQUIRE(baseArea.Footprint.Vertices[0].y == Approx(189.76460300106555));
    REQUIRE(baseArea.Buildings.size() == 1);
    REQUIRE(baseArea.Properties.size() == 2);

    Property property = baseArea.Properties[1];
    REQUIRE(property.FNR == 140127538);
    REQUIRE(property.UUID == "b8574e12-4618-4f1e-91fd-71e2f89b2375");
    REQUIRE(property.Buildings.empty());
    REQUIRE(property.Footprint.Vertices[0].x == Approx(539.78));
    REQUIRE(property.Footprint.Vertices[0].y == Approx(60.7778956));

    Building building = baseArea.Buildings[0];
    REQUIRE(building.UUID == "c8374e11-2767-4f0a-91fd-71d7f89b6681");
    REQUIRE(building.PropertyUUID == building.UUID);
    REQUIRE(building.PropertyFNR == 140029233);
    REQUIRE(building.BaseAreaID == "60605");
    REQUIRE(building.GroundHeight == Approx(34.865));
    REQUIRE(building.Height == Approx(107.574735807257));
    REQUIRE(building.error == 0);
    REQUIRE(building.Footprint.Vertices[0].x == Approx(551.020999965025));
    REQUIRE(building.Footprint.Vertices[0].y == Approx(57.5619951048866));

    if (firstRun)
    {
      json.clear();
      JSON::Serialize(district, json);
      JSON::Write(json, fileName2, 4);
    }
  }
}
