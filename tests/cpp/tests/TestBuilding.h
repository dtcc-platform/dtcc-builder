#include "BuildingProcessor.h"
#include "model/Building.h"

using namespace DTCC_BUILDER;

TEST_CASE("Point converage")
{
  Building building;
  building.footprint.vertices.push_back(Vector2D(0, 0));
  building.footprint.vertices.push_back(Vector2D(2, 0));
  building.footprint.vertices.push_back(Vector2D(2, 2));
  building.footprint.vertices.push_back(Vector2D(0, 2));
  building.footprint.vertices.push_back(Vector2D(0, 0));

  building.roof_points.push_back(Vector3D(0.51, 0.51, 0));
  building.roof_points.push_back(Vector3D(0.52, 0.52, 0));
  building.roof_points.push_back(Vector3D(0.53, 0.53, 0));

  building.roof_points.push_back(Vector3D(1.5, 0.5, 0));
  building.roof_points.push_back(Vector3D(1.5, 1.5, 0));

  REQUIRE(BuildingProcessor::point_coverage(building, 1.0) ==
          Approx(3.0 / 4.0));
  REQUIRE(BuildingProcessor::point_coverage(building, 0.5) ==
          Approx(3.0 / 16.0));
}
