#include "GeoJSON.h"
#include "JSON.h"
#include "RoadNetworkGenerator.h"
#include "SHP.h"

using namespace DTCC;

TEST_CASE("RoadNetwork generation")
{
  std::string filename1 = RootPath + "data/roadnetwork-torslanda/vl_riks.shp";
  std::string filename2 = RootPath + "data/roadnetwork-central-gbg/vl_riks.shp";

  std::vector<Polygon> roads;
  nlohmann::json attributes;

  // Test road data extracted from SHP.h
  SHP::Read(roads, filename1, nullptr, nullptr, &attributes);
  REQUIRE(roads.size() == 7);
  REQUIRE(roads[0].Vertices.size() == 13);
  Point2D v = roads[0].Vertices[0];
  REQUIRE(v.x == Approx(306234.751));
  REQUIRE(v.y == Approx(6401785.2819999997));

  json jsonAttrib = attributes["attributes"];
  REQUIRE(jsonAttrib.size() == 7);
  nlohmann::json firstAttr = jsonAttrib[0];
  REQUIRE(firstAttr.size() == 2);
  REQUIRE(firstAttr["KATEGORI"] == "Bilväg");
  REQUIRE(firstAttr["KKOD"] == "5071");

  // Check correct number of attributes in case of multi-part roads
  SHP::Read(roads, filename2, nullptr, nullptr, &attributes);
  REQUIRE(attributes["polyToAttribute"].size() == roads.size());

  // Test RoadNetwork object
  RoadNetwork roadNetwork = RoadNetworkGenerator::GetRoadNetwork(filename1);

  REQUIRE(roadNetwork.Edges.size() == 28);
  REQUIRE(roadNetwork.Vertices[0].x == v.x);
  REQUIRE(roadNetwork.Vertices[0].y == v.y);
  REQUIRE(roadNetwork.EdgeProperties.size() == firstAttr.size());
  REQUIRE(roadNetwork.EdgeProperties["KKOD"][0] == firstAttr["KKOD"]);
}

TEST_CASE("Convert GeoJSON to RoadNetwork")
{
  // Test writing to RoadNetwork.json
  // Smaller GeoJSON
  std::string path = RootPath + "data/geojson-to-roadnetwork/";
  GeoJSON::WriteRoadNetwork(path + "RoadNetwork.geojson");
  nlohmann::json jsonRoadNetwork;
  JSON::Read(jsonRoadNetwork, path + "RoadNetwork.json");

  REQUIRE(jsonRoadNetwork["RoadNetwork"]["Edges"].size() == 3);
  REQUIRE(jsonRoadNetwork["RoadNetwork"]["Vertices"].size() == 8);
  REQUIRE(jsonRoadNetwork["RoadNetwork"]["Vertices"][0] == 0.0);
  REQUIRE(jsonRoadNetwork["RoadNetwork"]["Vertices"][1] == Approx(1.9268));

  auto jsonValues = jsonRoadNetwork["RoadNetwork"]["EdgeValues"];
  REQUIRE(jsonValues["AI_w5km"][0] == Approx(0.89656978845596313));
  REQUIRE(jsonValues["AI_w2k_h"][2] == Approx(743.246826171875));

  // Larger GeoJSON
  GeoJSON::WriteRoadNetwork(path + "GOT_Network_Segmentmap.geojson");
  JSON::Read(jsonRoadNetwork, path + "GOT_Network_Segmentmap.json");

  json jsonNetwork = jsonRoadNetwork["RoadNetwork"];
  REQUIRE(jsonNetwork["Edges"].size() == 128349);
  REQUIRE(jsonNetwork["EdgeProperties"]["ID"][0] == "626044");
  REQUIRE(jsonNetwork["EdgeValues"]["AI_w2k_wl_h"][128348] ==
          Approx(2169.01025390625));

  // Test reading GeoJSON and creating RoadNetwork object
  RoadNetwork roadNetwork =
      GeoJSON::GetRoadNetwork(path + "RoadNetwork.geojson");
  REQUIRE(roadNetwork.Edges.size() == 3);
  REQUIRE(roadNetwork.Vertices[0].x == 101003.3176);
  REQUIRE(roadNetwork.EdgeValues["AI_w5km"][0] == Approx(0.89656978845596313));
  REQUIRE(roadNetwork.EdgeValues["AI_w2k_h"][2] == Approx(743.246826171875));
}