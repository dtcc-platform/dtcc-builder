#include "Polyfix.h"
#include "model/Building.h"
#include "model/City.h"
#include "model/Polygon.h"

using namespace DTCC_BUILDER;

TEST_CASE("Filter City")
{
  City base_model;
  Building building1;
  Building building2;
  Building building3;

  Polygon fp1;
  Polygon fp2;
  Polygon fp3;

  fp1.vertices.push_back(Vector2D(0, 0));
  fp1.vertices.push_back(Vector2D(0, 1));
  fp1.vertices.push_back(Vector2D(1, 1));
  fp1.vertices.push_back(Vector2D(1, 0));
  fp1.vertices.push_back(Vector2D(0, 0));

  Polyfix::make_closed(fp1, 0);
  Polyfix::make_oriented(fp1);

  fp2.vertices.push_back(Vector2D(0, 0));
  fp2.vertices.push_back(Vector2D(0, 5));
  fp2.vertices.push_back(Vector2D(5, 5));
  fp2.vertices.push_back(Vector2D(5, 0));
  fp2.vertices.push_back(Vector2D(0, 0));

  Polyfix::make_closed(fp2, 0);
  Polyfix::make_oriented(fp2);

  fp3.vertices.push_back(Vector2D(0, 0));
  fp3.vertices.push_back(Vector2D(0, 10));
  fp3.vertices.push_back(Vector2D(10, 10));
  fp3.vertices.push_back(Vector2D(10, 0));
  fp3.vertices.push_back(Vector2D(0, 0));

  Polyfix::make_closed(fp3, 0);
  Polyfix::make_oriented(fp3);

  building1.footprint = fp1;
  building1.error |= BuildingError::BUILDING_TOO_SMALL;
  building2.footprint = fp2;
  building2.error |= BuildingError::BUILDING_TOO_FEW_POINTS;
  building3.footprint = fp3;
  base_model.buildings.push_back(building1);
  base_model.buildings.push_back(building2);
  base_model.buildings.push_back(building3);

  /*
  City filteredModel;
  CityProcessor::BuildingFootprintFilter(base_model, filteredModel, 30);
  REQUIRE(base_model.buildings.size() == 3);
  REQUIRE(filteredModel.buildings.size() == 1);
  REQUIRE(Geometry::polygon_area(filteredModel.buildings[0].footprint) > 30);

  size_t tooSmallError = BuildingError::BUILDING_TOO_SMALL;
  size_t tooSmallFewPointsError = (BuildingError::BUILDING_TOO_SMALL |
                                   BuildingError::BUILDING_TOO_FEW_POINTS);
  size_t aspectError = BuildingError::BUILDING_BAD_ASPECT_RATIO;

  filteredModel.buildings.clear();
  CityProcessor::ErrorFilter(base_model, filteredModel, 0);
  REQUIRE(filteredModel.buildings.size() == 3);

  filteredModel.buildings.clear();
  CityProcessor::ErrorFilter(base_model, filteredModel, aspectError);
  REQUIRE(filteredModel.buildings.size() == 3);

  filteredModel.buildings.clear();
  CityProcessor::ErrorFilter(base_model, filteredModel, tooSmallError);
  REQUIRE(filteredModel.buildings.size() == 2);

  filteredModel.buildings.clear();
  CityProcessor::ErrorFilter(base_model, filteredModel,
                                  tooSmallFewPointsError);
  REQUIRE(filteredModel.buildings.size() == 1);
  */
}
