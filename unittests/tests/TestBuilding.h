#include "BuildingProcessor.h"
#include "datamodel/Building.h"

using namespace DTCC_BUILDER;

TEST_CASE("Point converage")
{
  Building building;
  building.Footprint.Vertices.push_back(Point2D(0, 0));
  building.Footprint.Vertices.push_back(Point2D(2, 0));
  building.Footprint.Vertices.push_back(Point2D(2, 2));
  building.Footprint.Vertices.push_back(Point2D(0, 2));
  building.Footprint.Vertices.push_back(Point2D(0, 0));

  building.RoofPoints.push_back(Point3D(0.51, 0.51, 0));
  building.RoofPoints.push_back(Point3D(0.52, 0.52, 0));
  building.RoofPoints.push_back(Point3D(0.53, 0.53, 0));

  building.RoofPoints.push_back(Point3D(1.5, 0.5, 0));
  building.RoofPoints.push_back(Point3D(1.5, 1.5, 0));

  REQUIRE(BuildingProcessor::PointCoverage(building, 1.0) == Approx(3.0 / 4.0));
  REQUIRE(BuildingProcessor::PointCoverage(building, 0.5) ==
          Approx(3.0 / 16.0));
}