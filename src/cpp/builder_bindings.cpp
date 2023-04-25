#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <pybind11/stl.h>

#include "CityModelGenerator.h"
#include "ElevationModelGenerator.h"
#include "GridField.h"
#include "LaplacianSmoother.h"
#include "Mesh.h"
#include "MeshGenerator.h"
#include "MeshProcessor.h"
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

GridField2D createBuilderGridField(py::array_t<double> data,
                                   py::tuple bounds,
                                   size_t XStep,
                                   size_t YStep,
                                   double XSize,
                                   double YSize)
{
  GridField2D gridField;
  double px = bounds[0].cast<double>();
  double py = bounds[1].cast<double>();
  double qx = bounds[2].cast<double>();
  double qy = bounds[3].cast<double>();
  auto bbox = BoundingBox2D(Point2D(px, py), Point2D(qx, qy));

  gridField.Grid.BoundingBox = bbox;
  gridField.Grid.XStep = XStep;
  gridField.Grid.YStep = YStep;
  gridField.Grid.XSize = XSize;
  gridField.Grid.YSize = YSize;

  auto data_r = data.unchecked<1>();
  size_t data_count = data_r.size();

  for (size_t i = 0; i < data_count; i++)
  {
    gridField.Values.push_back(data_r(i));
  }

  return gridField;
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

CityModel SimplifyCityModel(CityModel &cityModel,
                            py::tuple bounds,
                            double minimalBuildingDistance,
                            double minimalVertexDistance)
{
  double px = bounds[0].cast<double>();
  double py = bounds[1].cast<double>();
  double qx = bounds[2].cast<double>();
  double qy = bounds[3].cast<double>();
  auto bbox = BoundingBox2D(Point2D(px, py), Point2D(qx, qy));
  CityModelGenerator::SimplifyCityModel(
      cityModel, bbox, minimalBuildingDistance, minimalVertexDistance);
  return cityModel;
}

CityModel CleanCityModel(CityModel &cityModel, double minVertDistance)
{
  CityModelGenerator::CleanCityModel(cityModel, minVertDistance);
  return cityModel;
}

// Meshing

Mesh2D
GenerateMesh2D(const CityModel &cityModel, py::tuple bounds, double resolution)
{
  Mesh2D mesh;
  double px = bounds[0].cast<double>();
  double py = bounds[1].cast<double>();
  double qx = bounds[2].cast<double>();
  double qy = bounds[3].cast<double>();
  auto bbox = BoundingBox2D(Point2D(px, py), Point2D(qx, qy));

  MeshGenerator::GenerateMesh2D(mesh, cityModel, bbox, resolution);

  return mesh;
}

Mesh3D
GenerateMesh3D(const Mesh2D &mesh2D, double domainHeight, double meshResolution)
{
  Mesh3D mesh;
  auto num_layers =
      MeshGenerator::GenerateMesh3D(mesh, mesh2D, domainHeight, meshResolution);
  mesh.NumLayers = num_layers;
  return mesh;
}

// _pybuilder.SmoothMesh3D(mesh,city_model,dem,top_height,fix_buildings)
Mesh3D SmoothMesh3D(Mesh3D &mesh3D,
                    const CityModel &cityModel,
                    const GridField2D &dem,
                    bool fixBuildings)
{
  double topHeight{};
  LaplacianSmoother::SmoothMesh3D(mesh3D, cityModel, dem, topHeight,
                                  fixBuildings, false);
  return mesh3D;
}

std::vector<Surface3D> GenerateSurfaces3D(const CityModel &cityModel,
                                          const GridField2D &dtm,
                                          double resolution)
{
  Surface3D groundSurface;
  std::vector<Surface3D> buildingSurfaces;
  MeshGenerator::GenerateSurfaces3D(groundSurface, buildingSurfaces, cityModel,
                                    dtm, resolution);
  buildingSurfaces.insert(buildingSurfaces.begin(), groundSurface);
  return buildingSurfaces;
}

Mesh3D TrimMesh3D(Mesh3D &mesh3D,
                  const Mesh2D &mesh2D,
                  const CityModel &cityModel,
                  size_t numLayers)
{
  MeshGenerator::TrimMesh3D(mesh3D, mesh2D, cityModel, numLayers);
  return mesh3D;
}

Surface3D ExtractBoundary3D(const Mesh3D &mesh)
{
  Surface3D surface;
  MeshProcessor::ExtractBoundary3D(surface, mesh);
  return surface;
}

Surface3D ExtractOpenSurface3D(const Surface3D &boundary)
{
  Surface3D surface;
  MeshProcessor::ExtractOpenSurface3D(surface, boundary);
  return surface;
}

Surface3D MergeSurfaces3D(const std::vector<Surface3D> &surfaces)
{
  Surface3D merged_surface;
  MeshProcessor::MergeSurfaces3D(merged_surface, surfaces);
  return merged_surface;
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

  py::class_<DTCC_BUILDER::Vector3D>(m, "Vector3D")
      .def(py::init<>())
      .def_readonly("x", &DTCC_BUILDER::Vector3D::x)
      .def_readonly("y", &DTCC_BUILDER::Vector3D::y)
      .def_readonly("z", &DTCC_BUILDER::Vector3D::z);

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

  py::class_<DTCC_BUILDER::Simplex2D>(m, "Simplex2D")
      .def(py::init<>())
      .def_readonly("v0", &DTCC_BUILDER::Simplex2D::v0)
      .def_readonly("v1", &DTCC_BUILDER::Simplex2D::v1)
      .def_readonly("v2", &DTCC_BUILDER::Simplex2D::v2);

  py::class_<DTCC_BUILDER::Simplex3D>(m, "Simplex3D")
      .def(py::init<>())
      .def_readonly("v0", &DTCC_BUILDER::Simplex3D::v0)
      .def_readonly("v1", &DTCC_BUILDER::Simplex3D::v1)
      .def_readonly("v2", &DTCC_BUILDER::Simplex3D::v2)
      .def_readonly("v3", &DTCC_BUILDER::Simplex3D::v3);

  py::class_<DTCC_BUILDER::Mesh2D>(m, "Mesh2D").def(py::init<>());

  py::class_<DTCC_BUILDER::Mesh3D>(m, "Mesh3D")
      .def(py::init<>())
      .def_readonly("numLayers", &DTCC_BUILDER::Mesh3D::NumLayers)
      .def_readonly("Vertices", &DTCC_BUILDER::Mesh3D::Vertices)
      .def_readonly("Cells", &DTCC_BUILDER::Mesh3D::Cells)
      .def_readonly("Markers", &DTCC_BUILDER::Mesh3D::Markers);

  py::class_<DTCC_BUILDER::Surface3D>(m, "Surface3D")
      .def(py::init<>())
      .def_readonly("Vertices", &DTCC_BUILDER::Surface3D::Vertices)
      .def_readonly("Faces", &DTCC_BUILDER::Surface3D::Faces)
      .def_readonly("Normals", &DTCC_BUILDER::Surface3D::Normals);

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

  m.def("SimplifyCityModel", &DTCC_BUILDER::SimplifyCityModel,
        "Simplify city model");

  m.def("CleanCityModel", &DTCC_BUILDER::CleanCityModel, "Clean city model");

  m.def("GenerateMesh2D", &DTCC_BUILDER::GenerateMesh2D, "Generate 2D mesh");
  m.def("GenerateMesh3D", &DTCC_BUILDER::GenerateMesh3D, "Generate 2D mesh");

  m.def("SmoothMesh3D", &DTCC_BUILDER::SmoothMesh3D, "Smooth 3D mesh");
  m.def("TrimMesh3D", &DTCC_BUILDER::TrimMesh3D, "Trim 3D mesh");

  // m.def("GenerateSurface3D", &DTCC_BUILDER::GenerateSurface3D,
  //     "Generate 3D surface");
}
