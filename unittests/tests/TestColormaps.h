#include "Color.h"
#include "ColorMap.h"
#include "ColorMapIO.h"
#include "JSON.h"

using namespace DTCC;

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