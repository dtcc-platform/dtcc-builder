#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <pybind11/stl.h>

#include "CityModelGenerator.h"
#include "ElevationModelGenerator.h"
#include "GridField.h"
#include "Point.h"
#include "PointCloud.h"
#include "PointCloudMethods.h"
#include "Polygon.h"
#include "VertexSmoother.h"
#include "datamodel/Building.h"
#include "datamodel/CityModel.h"

namespace py = pybind11;

namespace DTCC_BUILDER
{
CityModel createBuilderCityModel(py::list footprints,
                                 py::list uuids,
                                 py::list heights,
                                 py::list ground_levels,
                                 py::tuple origin)
{
  CityModel cityModel;

  cityModel.Origin =
      Point2D(origin[0].cast<double>(), origin[1].cast<double>());
  size_t num_buildings = footprints.size();
  for (size_t i = 0; i < num_buildings; i++)
  {
    auto footprint = footprints[i].cast<py::list>();
    auto uuid = uuids[i].cast<std::string>();
    auto height = heights[i].cast<double>();
    auto ground_level = ground_levels[i].cast<double>();
    Polygon poly;
    for (size_t j = 0; j < footprint.size(); j++)
    {
      auto pt = footprint[j].cast<py::tuple>();
      poly.Vertices.push_back(
          Point2D(pt[0].cast<double>(), pt[1].cast<double>()));
    }
    Building building;
    building.Footprint = poly;
    building.UUID = uuid;
    building.Height = height;
    building.GroundHeight = ground_level;
    cityModel.Buildings.push_back(building);
  }
  return cityModel;
}

PointCloud createBuilderPointCloud(py::array_t<double> pts,
                                   py::array_t<uint8_t> cls,
                                   py::array_t<uint8_t> retNumber,
                                   py::array_t<uint8_t> numReturns)
{
  auto pts_r = pts.unchecked<2>();
  auto cls_r = cls.unchecked<1>();
  auto retNumber_r = retNumber.unchecked<1>();
  auto numReturns_r = numReturns.unchecked<1>();

  size_t pt_count = pts_r.shape(0);

  bool hasClassification = cls_r.size() == pt_count;
  bool hasReturnNumber = retNumber_r.size() == pt_count;
  bool hasNumberOfReturns = numReturns_r.size() == pt_count;

  PointCloud pointCloud;
  for (size_t i = 0; i < pt_count; i++)
  {
    pointCloud.Points.push_back(Point3D(pts_r(i, 0), pts_r(i, 1), pts_r(i, 2)));
    if (hasClassification)
      pointCloud.Classifications.push_back(cls_r(i));
    else
      pointCloud.Classifications.push_back(1);
    if (hasReturnNumber)
      pointCloud.ScanFlags.push_back(
          PointCloudProcessor::packScanFlag(retNumber_r(i), numReturns_r(i)));
  }

  pointCloud.BuildHasClassifications();
  pointCloud.CalculateBoundingBox();
  return pointCloud;
}

PointCloud removeGlobalOutliers(PointCloud &pointCloud, double outlierMargin)
{
  PointCloudProcessor::RemoveOutliers(pointCloud, outlierMargin);
  return pointCloud;
}

PointCloud removeVegetation(PointCloud &pointCloud)
{
  PointCloudProcessor::NaiveVegetationFilter(pointCloud);
  return pointCloud;
}

CityModel extractRoofPoints(CityModel &cityModel,
                            PointCloud &pointCloud,
                            double groundMargin,
                            double groundOutlierMargin,
                            double roofOutlierMargin,
                            size_t rootOutlerNeighbours,
                            double roofRANSACOutlierMargin,
                            size_t roofRANSACIterations)
{
  pointCloud = removeVegetation(pointCloud);
  CityModelGenerator::ExtractBuildingPoints(cityModel, pointCloud, groundMargin,
                                            groundOutlierMargin);
  if (roofOutlierMargin > 0)
  {
    CityModelGenerator::BuildingPointsOutlierRemover(
        cityModel, rootOutlerNeighbours, roofOutlierMargin);
  }
  if (roofRANSACIterations > 0)
  {
    CityModelGenerator::BuildingPointsRANSACOutlierRemover(
        cityModel, roofRANSACOutlierMargin, roofRANSACIterations);
  }
  return cityModel;
}

// GridField
GridField2D GenerateElevationModel(const PointCloud &pointCloud,
                                   double resolution,
                                   std::vector<int> classifications)
{
  GridField2D dem;
  ElevationModelGenerator::GenerateElevationModel(dem, pointCloud,
                                                  classifications, resolution);
  return dem;
}

GridField2D SmoothElevation(GridField2D &dem, size_t numSmoothings)
{
  VertexSmoother::SmoothField(dem, numSmoothings);
  return dem;
}

} // namespace DTCC_BUILDER

PYBIND11_MODULE(_pybuilder, m)
{
  py::class_<DTCC_BUILDER::CityModel>(m, "CityModel")
      .def(py::init<>())
      .def("__len__", [](const DTCC_BUILDER::CityModel &cm)
           { return cm.Buildings.size(); })
      .def_readonly("buildings", &DTCC_BUILDER::CityModel::Buildings)
      .def_readonly("origin", &DTCC_BUILDER::CityModel::Origin);

  py::class_<DTCC_BUILDER::Building>(m, "Building")
      .def(py::init<>())
      .def_readwrite("error", &DTCC_BUILDER::Building::error)
      .def_readwrite("uuid", &DTCC_BUILDER::Building::UUID)
      .def_readwrite("propertyUUID", &DTCC_BUILDER::Building::PropertyUUID)
      .def_readwrite("height", &DTCC_BUILDER::Building::Height)
      .def_readwrite("groundHeight", &DTCC_BUILDER::Building::GroundHeight)
      .def_readonly("footprint", &DTCC_BUILDER::Building::Footprint)
      .def_readonly("ground_points", &DTCC_BUILDER::Building::GroundPoints)
      .def_readonly("roof_points", &DTCC_BUILDER::Building::RoofPoints);

  py::class_<DTCC_BUILDER::Point2D>(m, "Point2D")
      .def(py::init<>())
      .def("__repr__",
           [](const DTCC_BUILDER::Point3D &p)
           {
             return "<Point3D (" + DTCC_BUILDER::str(p.x) + ", " +
                    DTCC_BUILDER::str(p.y) + ")>";
           })
      .def_readonly("x", &DTCC_BUILDER::Point2D::x)
      .def_readonly("y", &DTCC_BUILDER::Point2D::y);

  py::class_<DTCC_BUILDER::Point3D>(m, "Point3D")
      .def(py::init<>())
      .def("__repr__",
           [](const DTCC_BUILDER::Point3D &p)
           {
             return "<Point3D (" + DTCC_BUILDER::str(p.x) + ", " +
                    DTCC_BUILDER::str(p.y) + ", " + DTCC_BUILDER::str(p.z) +
                    ")>";
           })
      .def_readonly("x", &DTCC_BUILDER::Point3D::x)
      .def_readonly("y", &DTCC_BUILDER::Point3D::y)
      .def_readonly("z", &DTCC_BUILDER::Point3D::z);

  py::class_<DTCC_BUILDER::BoundingBox2D>(m, "BoundingBox")
      .def(py::init<>())
      .def_readonly("P", &DTCC_BUILDER::BoundingBox2D::P)
      .def_readonly("Q", &DTCC_BUILDER::BoundingBox2D::Q);

  py::class_<DTCC_BUILDER::Polygon>(m, "Polygon")
      .def(py::init<>())
      .def_readonly("vertices", &DTCC_BUILDER::Polygon::Vertices);

  py::class_<DTCC_BUILDER::PointCloud>(m, "PointCloud")
      .def(py::init<>())
      .def("__len__",
           [](const DTCC_BUILDER::PointCloud &p) { return p.Points.size(); })
      .def_readonly("points", &DTCC_BUILDER::PointCloud::Points)
      .def_readonly("classifications",
                    &DTCC_BUILDER::PointCloud::Classifications)
      .def_readonly("intensities", &DTCC_BUILDER::PointCloud::Intensities)
      .def_readonly("scan_flags", &DTCC_BUILDER::PointCloud::ScanFlags);

  py::class_<DTCC_BUILDER::GridField2D>(m, "GridField2D")
      .def(py::init<>())
      .def_readonly("Grid", &DTCC_BUILDER::GridField2D::Grid)
      .def_readonly("values", &DTCC_BUILDER::GridField2D::Values);

  py::class_<DTCC_BUILDER::Grid2D>(m, "Grid2D")
      .def(py::init<>())
      .def_readonly("xsize", &DTCC_BUILDER::Grid2D::XSize)
      .def_readonly("ysize", &DTCC_BUILDER::Grid2D::YSize)
      .def_readonly("xstep", &DTCC_BUILDER::Grid2D::XStep)
      .def_readonly("ystep", &DTCC_BUILDER::Grid2D::YStep);

  m.def("createBuilderCityModel", &DTCC_BUILDER::createBuilderCityModel,
        "create builder point cloud from citymodel data");

  m.def("createBuilderPointCloud", &DTCC_BUILDER::createBuilderPointCloud,
        "create builder point cloud from numpy arrays");
  m.def("removeGlobalOutliers", &DTCC_BUILDER::removeGlobalOutliers,
        "remove global outliers from point cloud");
  m.def("removeVegetation", &DTCC_BUILDER::removeVegetation,
        "remove vegetation from point cloud");
  m.def("extractRoofPoints", &DTCC_BUILDER::extractRoofPoints,
        "extract roof points from point cloud");
  m.def("GenerateElevationModel", &DTCC_BUILDER::GenerateElevationModel,
        "generate height field from point cloud");
  m.def("SmoothElevation", &DTCC_BUILDER::SmoothElevation,
        "Smooth  elevation grid field");
}
