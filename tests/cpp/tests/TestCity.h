#include "Polyfix.h"
#include "model/Building.h"
#include "model/City.h"
#include "model/Polygon.h"

using namespace DTCC_BUILDER;

TEST_CASE("Filter City")
{
  City baseModel;
  Building building1;
  Building building2;
  Building building3;

  Polygon fp1;
  Polygon fp2;
  Polygon fp3;

  fp1.Vertices.push_back(Point2D(0, 0));
  fp1.Vertices.push_back(Point2D(0, 1));
  fp1.Vertices.push_back(Point2D(1, 1));
  fp1.Vertices.push_back(Point2D(1, 0));
  fp1.Vertices.push_back(Point2D(0, 0));

  Polyfix::MakeClosed(fp1, 0);
  Polyfix::MakeOriented(fp1);

  fp2.Vertices.push_back(Point2D(0, 0));
  fp2.Vertices.push_back(Point2D(0, 5));
  fp2.Vertices.push_back(Point2D(5, 5));
  fp2.Vertices.push_back(Point2D(5, 0));
  fp2.Vertices.push_back(Point2D(0, 0));

  Polyfix::MakeClosed(fp2, 0);
  Polyfix::MakeOriented(fp2);

  fp3.Vertices.push_back(Point2D(0, 0));
  fp3.Vertices.push_back(Point2D(0, 10));
  fp3.Vertices.push_back(Point2D(10, 10));
  fp3.Vertices.push_back(Point2D(10, 0));
  fp3.Vertices.push_back(Point2D(0, 0));

  Polyfix::MakeClosed(fp3, 0);
  Polyfix::MakeOriented(fp3);

  building1.Footprint = fp1;
  building1.error |= BuildingError::BUILDING_TOO_SMALL;
  building2.Footprint = fp2;
  building2.error |= BuildingError::BUILDING_TOO_FEW_POINTS;
  building3.Footprint = fp3;
  baseModel.Buildings.push_back(building1);
  baseModel.Buildings.push_back(building2);
  baseModel.Buildings.push_back(building3);

  /*
  City filteredModel;
  CityProcessor::BuildingFootprintFilter(baseModel, filteredModel, 30);
  REQUIRE(baseModel.Buildings.size() == 3);
  REQUIRE(filteredModel.Buildings.size() == 1);
  REQUIRE(Geometry::PolygonArea(filteredModel.Buildings[0].Footprint) > 30);

  size_t tooSmallError = BuildingError::BUILDING_TOO_SMALL;
  size_t tooSmallFewPointsError = (BuildingError::BUILDING_TOO_SMALL |
                                   BuildingError::BUILDING_TOO_FEW_POINTS);
  size_t aspectError = BuildingError::BUILDING_BAD_ASPECT_RATIO;

  filteredModel.Buildings.clear();
  CityProcessor::ErrorFilter(baseModel, filteredModel, 0);
  REQUIRE(filteredModel.Buildings.size() == 3);

  filteredModel.Buildings.clear();
  CityProcessor::ErrorFilter(baseModel, filteredModel, aspectError);
  REQUIRE(filteredModel.Buildings.size() == 3);

  filteredModel.Buildings.clear();
  CityProcessor::ErrorFilter(baseModel, filteredModel, tooSmallError);
  REQUIRE(filteredModel.Buildings.size() == 2);

  filteredModel.Buildings.clear();
  CityProcessor::ErrorFilter(baseModel, filteredModel,
                                  tooSmallFewPointsError);
  REQUIRE(filteredModel.Buildings.size() == 1);
  */
}
