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
#include "Mesh.h"
#include "MeshField.h"
#include "MeshVectorField.h"


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

TEST_CASE("COLORMAPS")
{
  ColorMap cm;
  cm.InsertColor(1,Color(1.0,1.0,1.0));
  cm.InsertColor(0,Color(0.0,0.0,0.0));

  ColorMap cm2;
  cm2.InsertColor(0.8,Color(0.0,0.0,0.0));
  cm2.InsertColor(0.9,Color(1.0,0.0,0.0));

  SECTION("Insert")
  {
    REQUIRE(cm.size() == 2);
    REQUIRE(cm.Colors.front().first == 0);
    REQUIRE(cm.Colors.back().first == 1);
  }
  SECTION("Interpolate")
  {
    REQUIRE(cm(1.1).R == 1.0);
    REQUIRE(cm(1.1).G == 1.0);
    REQUIRE(cm(1.1).B == 1.0);

    REQUIRE(cm(0.5).R == 0.5);
    REQUIRE(cm(0.5).G == 0.5);
    REQUIRE(cm(0.5).B == 0.5);

    REQUIRE(cm2(0.85).R == Approx(0.5).margin(0.0001));
    REQUIRE(cm2(0.85).G == 0.0);
    REQUIRE(cm2(0.85).B == 0.0);
  }

  SECTION("Load PNG")
  { 
    ColorMap cm3;
    ColorMapIO::ReadPNG(cm3, "../unittests/data/colormap_jet.png");
    REQUIRE(cm3.size()==256);
    REQUIRE(cm3(0).R==Approx(0).margin(0.0001));
    REQUIRE(cm3(0).G==Approx(0).margin(0.0001));
    REQUIRE(cm3(0).B==Approx(127/255.0).margin(0.0001));

    REQUIRE(cm3(0.5).R==Approx(121/255.0).margin(0.0001));
    REQUIRE(cm3(0.5).G==Approx(255/255.0).margin(0.0001));
    REQUIRE(cm3(0.5).B==Approx(124.5/255.0).margin(0.0001));

    REQUIRE(cm3(1).R==Approx(127/255.0).margin(0.0001));
    REQUIRE(cm3(1).G==Approx(0).margin(0.0001));
    REQUIRE(cm3(1).B==Approx(0).margin(0.0001));

  }
}
