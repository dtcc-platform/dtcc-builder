// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#define CATCH_CONFIG_MAIN

#include "Grid.h"
#include "GridField.h"
#include "GridVectorField.h"
#include "Mesh.h"
#include "MeshField.h"
#include "MeshVectorField.h"
#include "catch.hpp"
#include <XMLParser.h>
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
