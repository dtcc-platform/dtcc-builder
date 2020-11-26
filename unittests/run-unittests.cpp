// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#define CATCH_CONFIG_MAIN

#include "catch.hpp"

#include "Color.h"
#include "ColorMap.h"
#include "ColorMapIO.h"
#include "Grid.h"
#include "GridField.h"
#include "GridVectorField.h"
#include "Hashing.h"
#include "JSON.h"
#include "Mesh.h"
#include "MeshField.h"
#include "MeshVectorField.h"
#include <XMLParser.h>
#include <deque>
#include <nlohmann/json.hpp>

using namespace DTCC;

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
  const std::string RootPath1 =
      "/home/dtcc/core/unittests/xml-parser-test-files/";
  const std::string RootPath2 =
      "/builds/dtcc3d/core/unittests/xml-parser-test-files/";
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
    ColorMapIO::ReadPNG(cm3, "../unittests/data/colormap_jet.png");
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
    ColorMapIO::ReadCPT(cm6, "../unittests/data/inferno.cpt");
    REQUIRE(cm6.size() == 255 * 2);
    REQUIRE(cm6(125 / 255.0).R == Approx(183 / 255.0).margin(0.0001));
    REQUIRE(cm6(125 / 255.0).G == Approx(53 / 255.0).margin(0.0001));
    REQUIRE(cm6(125 / 255.0).B == Approx(87 / 255.0).margin(0.0001));

    ColorMap cm7;
    ColorMapIO::ReadCPT(cm7, "../unittests/data/BrBG_11.cpt");
    REQUIRE(cm7(0.5).R == Approx(245 / 255.0).margin(0.0001));
    REQUIRE(cm7(0.5).G == Approx(245 / 255.0).margin(0.0001));
    REQUIRE(cm7(0.5).B == Approx(245 / 255.0).margin(0.0001));
  }

  SECTION("Serialize JSON")
  {
    ColorMap cm3;
    ColorMapIO::ReadPNG(cm3, "../unittests/data/colormap_jet.png");
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

TEST_CASE("citymodel")
{
  const char *fileName1 = "../unittests/data/CityModelExample.json";
  const char *fileName2 = "../unittests/data/CityModelExample2.json";

  for (bool firstRun : {true, false})
  {
    District district;
    nlohmann::json json;

    JSON::Read(json, firstRun ? fileName1 : fileName2);
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
    REQUIRE(primaryArea.Footprint.Vertices[0].x == Approx(-446.4952344278572));
    REQUIRE(primaryArea.Footprint.Vertices[0].y == Approx(150.96198354940861));

    BaseArea baseArea = primaryArea.BaseAreas[0];
    REQUIRE(baseArea.AreaID == "60605");
    REQUIRE(baseArea.PrimaryAreaID == "606");
    REQUIRE(baseArea.Footprint.Vertices[0].x == Approx(214.59035302355187));
    REQUIRE(baseArea.Footprint.Vertices[0].y == Approx(189.76460300106555));
    REQUIRE(baseArea.Buildings.size() == 1);
    REQUIRE(baseArea.Properties.size() == 1);

    Property property = baseArea.Properties[0];
    REQUIRE(property.FNR == 140029233);
    REQUIRE(property.UUID == "c8374e11-2767-4f0a-91fd-71d7f89b6681");
    REQUIRE(property.Buildings.size() == 1);
    REQUIRE(property.Footprint.Vertices[0].x == Approx(529.5));
    REQUIRE(property.Footprint.Vertices[0].y == Approx(59.0230000000447));

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
