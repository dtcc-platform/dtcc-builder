// Copyright (C) 2020-2021 Anders Logg, Anton J Olsson
// Licensed under the MIT License

#define CATCH_CONFIG_MAIN

#include "CSV.h"
#include "Color.h"
#include "ColorMap.h"
#include "ColorMapIO.h"
#include "GeoJSON.h"
#include "Grid.h"
#include "GridField.h"
#include "GridVectorField.h"
#include "Hashing.h"
#include "JSON.h"
#include "LAS.h"
#include "Mesh.h"
#include "MeshField.h"
#include "MeshVectorField.h"
#include "PointCloud.h"
#include "PointCloudProcessor.h"
#include "SHP.h"
#include "catch.hpp"
#include <RoadNetworkGenerator.h>
#include <XMLParser.h>
#include <deque>
#include <nlohmann/json.hpp>

using namespace DTCC;

const std::string RootPath{"/home/dtcc/core/unittests/"};

TEST_CASE("Grid2D")
{
  Point2D p(0, 0);
  Point2D q(1, 1);
  BoundingBox2D bbox(p, q);
  Grid2D grid(bbox, 4, 5);

  SECTION("StepSize")
  {
    REQUIRE(grid.XStep == Approx(1.0/3.0));
    REQUIRE(grid.YStep == Approx(1.0/4.0));
  }

  SECTION("NumVertices")
  {
    REQUIRE(grid.NumVertices() == 20);
  }

  SECTION("NumCells")
  {
    REQUIRE(grid.NumCells() == 12);
  }

  SECTION("Index2Point2Index")
  {
    size_t index = grid.NumVertices() / 3;
    REQUIRE(grid.Point2Index(grid.Index2Point(index)) == index);
  }
}

TEST_CASE("Grid3D")
{
  Point3D p(0, 0, 0);
  Point3D q(1, 1, 1);
  BoundingBox3D bbox(p, q);
  Grid3D grid(bbox, 4, 5, 6);

  SECTION("StepSize")
  {
    REQUIRE(grid.XStep == Approx(1.0/3.0));
    REQUIRE(grid.YStep == Approx(1.0/4.0));
    REQUIRE(grid.ZStep == Approx(1.0/5.0));
  }

  SECTION("NumVertices")
  {
    REQUIRE(grid.NumVertices() == 120);
  }

  SECTION("NumCells")
  {
    REQUIRE(grid.NumCells() == 60);
  }

  SECTION("Index2Point2Index")
  {
    size_t index = grid.NumVertices() / 3;
    REQUIRE(grid.Point2Index(grid.Index2Point(index)) == index);
  }
}

TEST_CASE("GridField2D")
{
  Point2D p(0, 0);
  Point2D q(1, 1);
  Point2D r(0.234, 0.485);
  BoundingBox2D bbox(p, q);
  Grid2D grid(bbox, 11, 17);
  GridField2D u(grid);

  class MyField : public Field2D
  {
  public:
    double operator()(const Point2D& p) const
    {
      return sin(p.x)*cos(p.y);
    }
  };

  SECTION("Evaluate")
  {
    REQUIRE(u(p) == Approx(0.0));
    REQUIRE(u(q) == Approx(0.0));
  }

  SECTION("Interpolate")
  {
    MyField f;
    u.Interpolate(f);
    REQUIRE(u(r) == Approx(f(r)).margin(0.01));
  }
}

TEST_CASE("GridField3D")
{
  Point3D p(0, 0, 0);
  Point3D q(1, 1, 1);
  Point3D r(0.234, 0.485, 0.763);
  BoundingBox3D bbox(p, q);
  Grid3D grid(bbox, 11, 17, 23);
  GridField3D u(grid);

  class MyField : public Field3D
  {
  public:
    double operator()(const Point3D& p) const
    {
      return sin(p.x)*cos(p.y)*exp(p.z);
    }
  };

  SECTION("Evaluate")
  {
    REQUIRE(u(p) == Approx(0.0));
    REQUIRE(u(q) == Approx(0.0));
  }

  SECTION("Interpolate")
  {
    MyField f;
    u.Interpolate(f);
    REQUIRE(u(r) == Approx(f(r)).margin(0.01));
  }
}

TEST_CASE("GridVectorField2D")
{
  Point2D p(0, 0);
  Point2D q(1, 1);
  Point2D r(0.234, 0.485);
  BoundingBox2D bbox(p, q);
  Grid2D grid(bbox, 11, 17);
  GridVectorField2D u(grid);

  class MyField : public VectorField2D
  {
  public:
    Vector2D operator()(const Point2D& p) const
    {
      return Vector2D(sin(p.x)*cos(p.y), cos(p.x)*sin(p.y));
    }
  };

  SECTION("Evaluate")
  {
    REQUIRE(u(p).x == Approx(0.0));
    REQUIRE(u(p).y == Approx(0.0));
    REQUIRE(u(q).x == Approx(0.0));
    REQUIRE(u(q).y == Approx(0.0));
  }

  SECTION("Interpolate")
  {
    MyField f;
    u.Interpolate(f);
    REQUIRE(u(r).x == Approx(f(r).x).margin(0.01));
    REQUIRE(u(r).y == Approx(f(r).y).margin(0.01));
  }
}

TEST_CASE("GridVectorField3D")
{
  Point3D p(0, 0, 0);
  Point3D q(1, 1, 1);
  Point3D r(0.234, 0.485, 0.763);
  BoundingBox3D bbox(p, q);
  Grid3D grid(bbox, 11, 17, 23);
  GridVectorField3D u(grid);

  class MyField : public VectorField3D
  {
  public:
    Vector3D operator()(const Point3D& p) const
    {
      return Vector3D(sin(p.x)*cos(p.y)*exp(p.z),
                      cos(p.x)*exp(p.y)*sin(p.z),
                      exp(p.x)*sin(p.y)*cos(p.z));
    }
  };

  SECTION("Evaluate")
  {
    REQUIRE(u(p).x == Approx(0.0));
    REQUIRE(u(p).y == Approx(0.0));
    REQUIRE(u(p).z == Approx(0.0));
    REQUIRE(u(q).x == Approx(0.0));
    REQUIRE(u(q).y == Approx(0.0));
    REQUIRE(u(q).z == Approx(0.0));
  }

  SECTION("Interpolate")
  {
    MyField f;
    u.Interpolate(f);
    REQUIRE(u(r).x == Approx(f(r).x).margin(0.01));
    REQUIRE(u(r).y == Approx(f(r).y).margin(0.01));
    REQUIRE(u(r).z == Approx(f(r).z).margin(0.01));
  }
}

TEST_CASE("XMLParser")
{
  const std::string RootPath1 = RootPath + "data/xml-parser-test-files/";
  const std::string RootPath2 = RootPath + "data/xml-parser-test-files/";
  std::string rootPath = RootPath1;

  pugi::xml_document doc;
  std::string filePath = rootPath + "XMLExampleFile1.xml";
  pugi::xml_parse_result result = doc.load_file(filePath.c_str());
  if (!result)
    rootPath = RootPath2;
  filePath = rootPath + "XMLExampleFile1.xml";

  nlohmann::json json = XMLParser::GetJsonFromXML(filePath.c_str(), true);
  REQUIRE(json["city"]["name"] == "Johanneberg");

  json = XMLParser::GetJsonFromXML(filePath.c_str());
  REQUIRE(json["geometry"].is_array());
  REQUIRE(json["geometry"][0]["file_name"] == "Mesh/Buildings.vtk");

  json = XMLParser::GetJsonFromXML((rootPath + "XMLExampleFile2.xml").c_str());
  REQUIRE(json["Test"][0]["TestId"] == "0001");
  REQUIRE(json["example1"]["#content"][0] == 42);
  REQUIRE(json["example1"]["tag"][0] == "Content1");
  REQUIRE(json["example1"].size() == 2);
  REQUIRE(json["example1"]["#content"].size() == 2);
  REQUIRE(json["example2"]["#content"] == 37);
  REQUIRE(json["example2"]["x"] == 5);
  REQUIRE(json["example3"] == -26.3);

  json = XMLParser::GetJsonFromXML((rootPath + "XMLExampleFile3.xml").c_str(),
                                   false);
  REQUIRE(json["item"]["batters"]["batter"][0]["id"] == 1001);
  REQUIRE(json["item"]["topping"][5]["#content"] == "Maple");
}

TEST_CASE("COLORMAPS")
{
  ColorMap cm;
  cm.InsertColor(1, Color(1.0, 1.0, 1.0));
  cm.InsertColor(0, Color(0.0, 0.0, 0.0));

  ColorMap cm2;
  cm2.InsertColor(0.8, Color(0.0, 0.0, 0.0));
  cm2.InsertColor(0.9, Color(1.0, 0.0, 0.0));

  SECTION("Insert")
  {
    REQUIRE(cm.size() == 2);
    REQUIRE(cm.Colors.front().first == 0);
    REQUIRE(cm.Colors.back().first == 1);
  }
  SECTION("Interpolate")
  {

    REQUIRE(cm(0).R == 0.0);
    REQUIRE(cm(0).G == 0.0);
    REQUIRE(cm(0).B == 0.0);

    REQUIRE(cm(1.0).R == 1.0);
    REQUIRE(cm(1.0).G == 1.0);
    REQUIRE(cm(1.0).B == 1.0);

    REQUIRE(cm(1.1).R == 1.0);
    REQUIRE(cm(1.1).G == 1.0);
    REQUIRE(cm(1.1).B == 1.0);

    REQUIRE(cm(0.3).R == 0.3);
    REQUIRE(cm(0.3).G == 0.3);
    REQUIRE(cm(0.3).B == 0.3);

    REQUIRE(cm2(0.85).R == Approx(0.5).margin(0.0001));
    REQUIRE(cm2(0.85).G == 0.0);
    REQUIRE(cm2(0.85).B == 0.0);
  }

  SECTION("Load PNG")
  {
    ColorMap cm3;
    ColorMapIO::ReadPNG(cm3, RootPath + "data/colormap_jet.png");
    REQUIRE(cm3.size() == 256);

    REQUIRE(cm3(0).R == Approx(127 / 255.0).margin(0.0001));
    REQUIRE(cm3(0).G == Approx(0).margin(0.0001));
    REQUIRE(cm3(0).B == Approx(0).margin(0.0001));

    REQUIRE(cm3(0.5).R == Approx(121 / 255.0).margin(0.0001));
    REQUIRE(cm3(0.5).G == Approx(255 / 255.0).margin(0.0001));
    REQUIRE(cm3(0.5).B == Approx(124.5 / 255.0).margin(0.0001));

    REQUIRE(cm3(1).R == Approx(0).margin(0.0001));
    REQUIRE(cm3(1).G == Approx(0).margin(0.0001));
    REQUIRE(cm3(1).B == Approx(127 / 255.0).margin(0.0001));
  }

  SECTION("Write PNG")
  {
    ColorMapIO::WritePNG(cm, "testmap.png");
    ColorMap cm4;
    ColorMapIO::ReadPNG(cm4, "testmap.png");
    REQUIRE(cm4(0.3).R == Approx(0.3).margin(0.0001));
    REQUIRE(cm4(0.3).G == Approx(0.3).margin(0.0001));
    REQUIRE(cm4(0.3).B == Approx(0.3).margin(0.0001));
    remove("testmap.png");
  }

  SECTION("Read cpt")
  {
    ColorMap cm6;
    ColorMapIO::ReadCPT(cm6, RootPath + "data/inferno.cpt");
    REQUIRE(cm6.size() == 255 * 2);
    REQUIRE(cm6(125 / 255.0).R == Approx(183 / 255.0).margin(0.0001));
    REQUIRE(cm6(125 / 255.0).G == Approx(53 / 255.0).margin(0.0001));
    REQUIRE(cm6(125 / 255.0).B == Approx(87 / 255.0).margin(0.0001));

    ColorMap cm7;
    ColorMapIO::ReadCPT(cm7, RootPath + "data/BrBG_11.cpt");
    REQUIRE(cm7(0.5).R == Approx(245 / 255.0).margin(0.0001));
    REQUIRE(cm7(0.5).G == Approx(245 / 255.0).margin(0.0001));
    REQUIRE(cm7(0.5).B == Approx(245 / 255.0).margin(0.0001));
  }

  SECTION("Serialize JSON")
  {
    ColorMap cm3;
    ColorMapIO::ReadPNG(cm3, RootPath + "data/colormap_jet.png");
    JSON::Write(cm3, "testmap.json");
    ColorMap cm5;
    JSON::Read(cm5, "testmap.json");

    REQUIRE(cm5.size() == 256);

    REQUIRE(cm5(0).R == Approx(127 / 255.0).margin(0.0001));
    REQUIRE(cm5(0).G == Approx(0).margin(0.0001));
    REQUIRE(cm5(0).B == Approx(0).margin(0.0001));

    REQUIRE(cm5(0.5).R == Approx(121 / 255.0).margin(0.0001));
    REQUIRE(cm5(0.5).G == Approx(255 / 255.0).margin(0.0001));
    REQUIRE(cm5(0.5).B == Approx(124.5 / 255.0).margin(0.0001));

    REQUIRE(cm5(1).R == Approx(0).margin(0.0001));
    REQUIRE(cm5(1).G == Approx(0).margin(0.0001));
    REQUIRE(cm5(1).B == Approx(127 / 255.0).margin(0.0001));
    // remove("testmap.json");
  }
}

TEST_CASE("datamodel")
{
  const std::string fileName1 = RootPath + "data/CityModelExample.json";
  const std::string fileName2 = RootPath + "data/CityModelExample2.json";

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

TEST_CASE("Hashing")
{
  SECTION("Hash Point2D")
  {
    Point2D p(1, 2);
    Info(Hashing::Hex(Hashing::Hash(p)));
  }

  SECTION("Hash Point3D")
  {
    Point3D p(1, 2, 3);
    Info(Hashing::Hex(Hashing::Hash(p)));
  }
}

TEST_CASE("ISO 8559-1 to UTF-8")
{
  std::string testStr("G\345ngv\344g");
  REQUIRE(Utils::Iso88591ToUtf8(testStr) == "Gångväg");
}

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

TEST_CASE("POINT_CLOUD")
{
  SECTION("READ LAS")
  {
    PointCloud pc;
    LAS::Read(pc, RootPath + "data/minimal_las.las");

    REQUIRE(pc.Points.size()==10);
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
    pc.Points.push_back(Vector3D(0,0,0));
    pc.Classifications.push_back(0);
    pc.Points.push_back(Vector3D(1,0,0));
    pc.Classifications.push_back(1);
    pc.Points.push_back(Vector3D(2,0,0));
    pc.Classifications.push_back(2);

    PointCloud out_pc = PointCloudProcessor::ClassificationFilter(pc,{1,2});
    REQUIRE(out_pc.Points.size()==2);
    REQUIRE(out_pc.Classifications.size() == 2);
    REQUIRE(out_pc.Points[0].x==1);
    REQUIRE(out_pc.Points[1].x==2);

  }
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

  Point2D p1 = Point2D(0,0);
  Point2D p2 = Point2D(5,5);
  Point2D p3 = Point2D(-5,-5);
  Point2D p4 = Point2D(-1,-1);
  Point2D p5 = Point2D(-2,-2);
  Point2D p6 = Point2D(2,2);
  Point2D p7 = Point2D(1,1);

  SECTION("UNION")
  {
    BoundingBox2D bb1 = BoundingBox2D(p1,p2);
    BoundingBox2D bb2 = BoundingBox2D(p3,p4);
    bb1.Union(bb2);
    REQUIRE(bb1.P.x == -5);
    REQUIRE(bb1.P.y == -5);
    REQUIRE(bb1.Q.x == 5);
    REQUIRE(bb1.Q.y == 5);
  }

  SECTION("INTERSECT")
  {
    BoundingBox2D bb1 = BoundingBox2D(p1,p2);
    BoundingBox2D bb2 = BoundingBox2D(p3,p4);
    bb1.Intersect(bb2);
    REQUIRE(bb1.Area() == 0);

    BoundingBox2D bb3 = BoundingBox2D(p1,p2);
    BoundingBox2D bb4 = BoundingBox2D(p5,p6);
    bb3.Intersect(bb4);
    REQUIRE(bb3.P.x == 0);
    REQUIRE(bb3.P.y == 0);
    REQUIRE(bb3.Q.x == 2);
    REQUIRE(bb3.Q.y == 2);

    BoundingBox2D bb5= BoundingBox2D(p1,p2);
    BoundingBox2D bb6 = BoundingBox2D(p7,p6);
    bb5.Intersect(bb6);
    REQUIRE(bb5.P.x == bb6.P.x);
    REQUIRE(bb5.P.y == bb6.P.y);
    REQUIRE(bb5.Q.x == bb6.Q.x);
    REQUIRE(bb5.Q.y == bb6.Q.y);
  }

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
    REQUIRE(
        std::count_if(UUIDs.begin(), UUIDs.end(), [](const std::string &UUID) {
          return UUID.length() == defUUIDLength;
        }) == UUIDs.size());
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

TEST_CASE("Subtract polygon origin")
{
  std::vector<Point2D> vertices = {Point2D(3, 7), Point2D(2, 4.5)};
  Polygon p;
  p.Vertices = vertices;
  Polyfix::Transform(p, Point2D(1, 2));
  Polyfix::Transform(vertices, Point2D(1, 2));
  for (const auto &vertices2 : {vertices, p.Vertices})
  {
    REQUIRE(vertices2[0].x == 2);
    REQUIRE(vertices2[0].y == 5);
    REQUIRE(vertices2[1].x == 1);
    REQUIRE(vertices2[1].y == 2.5);
  }
}

TEST_CASE("Parse Parameters")
{
  Parameters p;
  p["DataDirectory"] = "359a538b-4616-4c33-b9b0-3870776fa28a";
  JSON::Write(p, "Parameters.json");
  Parameters q;
  JSON::Read(q, "Parameters.json");
  const std::string d = q["DataDirectory"];
  REQUIRE(d == "359a538b-4616-4c33-b9b0-3870776fa28a");
}

TEST_CASE("Get Filename from path")
{
  std::string path = "/path/to/fileName.json";
  std::string pathFileOnly = "fileName.json";
  std::string pathNoFile = "/path/to/dir/";

  REQUIRE(Utils::GetFilename(path) == "fileName.json");
  REQUIRE(Utils::GetFilename(path, true) == "fileName");
  REQUIRE(Utils::GetFilename(pathFileOnly) == "fileName.json");
  REQUIRE(Utils::GetFilename(pathFileOnly, true) == "fileName");
  REQUIRE(Utils::GetFilename(pathNoFile) == "");
}
