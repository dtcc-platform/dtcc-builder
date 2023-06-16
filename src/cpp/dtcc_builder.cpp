// Copyright (C) 2023 Dag Wästberg
// Licensed under the MIT License
//
// Modified by Anders Logg 2023

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <pybind11/stl.h>

#include "CityModelGenerator.h"
#include "ElevationModelGenerator.h"
#include "MeshGenerator.h"
#include "MeshProcessor.h"
#include "Smoother.h"
#include "VertexSmoother.h"
#include "model/Building.h"
#include "model/CityModel.h"
#include "model/GridField.h"
#include "model/Mesh.h"
#include "model/Point.h"
#include "model/PointCloud.h"
#include "model/Polygon.h"

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
  CityModelGenerator::CleanCityModel(cityModel, 1.0);

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

GridField createBuilderGridField(py::array_t<double> data,
                                 py::tuple bounds,
                                 size_t XSize,
                                 size_t YSize,
                                 double XStep,
                                 double YStep)
{
  GridField gridField;
  double px = bounds[0].cast<double>();
  double py = bounds[1].cast<double>();
  double qx = bounds[2].cast<double>();
  double qy = bounds[3].cast<double>();
  auto bbox = BoundingBox2D(Point2D(px, py), Point2D(qx, qy));

  gridField.grid.BoundingBox = bbox;
  gridField.grid.XStep = XStep;
  gridField.grid.YStep = YStep;
  gridField.grid.XSize = XSize;
  gridField.grid.YSize = YSize;

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
GridField GenerateElevationModel(const PointCloud &pointCloud,
                                 double resolution,
                                 std::vector<int> classifications)
{
  GridField dem;
  ElevationModelGenerator::GenerateElevationModel(dem, pointCloud,
                                                  classifications, resolution);
  return dem;
}

GridField SmoothElevation(GridField &dem, size_t numSmoothings)
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
      cityModel, bbox, minimalBuildingDistance, minimalVertexDistance / 2);
  return cityModel;
}

CityModel CleanCityModel(CityModel &cityModel, double minVertDistance)
{
  CityModelGenerator::CleanCityModel(cityModel, minVertDistance / 2);
  return cityModel;
}

// Meshing

Mesh GenerateMesh2D(const CityModel &cityModel,
                    py::tuple bounds,
                    double resolution)
{
  Mesh mesh;
  double px = bounds[0].cast<double>();
  double py = bounds[1].cast<double>();
  double qx = bounds[2].cast<double>();
  double qy = bounds[3].cast<double>();
  auto bbox = BoundingBox2D(Point2D(px, py), Point2D(qx, qy));

  MeshGenerator::GenerateMesh2D(mesh, cityModel, bbox, resolution);

  return mesh;
}

VolumeMesh
GenerateVolumeMesh(const Mesh &mesh, double domainHeight, double meshResolution)
{
  VolumeMesh volume_mesh;
  auto num_layers = MeshGenerator::GenerateVolumeMesh(
      volume_mesh, mesh, domainHeight, meshResolution);
  volume_mesh.NumLayers = num_layers;
  return volume_mesh;
}

VolumeMesh smooth_volume_mesh(VolumeMesh &volume_mesh,
                              const CityModel &city_model,
                              const GridField &dem,
                              double top_height,
                              bool fix_buildings,
                              size_t max_iterations,
                              double relative_tolerance)
{
  Smoother::smooth_volume_mesh(volume_mesh, city_model, dem, top_height,
                               fix_buildings, max_iterations,
                               relative_tolerance);
  return volume_mesh;
}

std::vector<Mesh> GenerateSurfaces3D(const CityModel &cityModel,
                                     const GridField &dtm,
                                     double resolution)
{
  Mesh groundSurface;
  std::vector<Mesh> buildingSurfaces;
  MeshGenerator::GenerateSurfaces3D(groundSurface, buildingSurfaces, cityModel,
                                    dtm, resolution);
  buildingSurfaces.insert(buildingSurfaces.begin(), groundSurface);
  return buildingSurfaces;
}

VolumeMesh TrimVolumeMesh(VolumeMesh &volume_mesh,
                          const Mesh &mesh,
                          const CityModel &cityModel)
{
  size_t numLayers = volume_mesh.NumLayers;
  MeshGenerator::TrimVolumeMesh(volume_mesh, mesh, cityModel, numLayers);
  return volume_mesh;
}

Mesh ExtractBoundary3D(const VolumeMesh &mesh)
{
  Mesh surface;
  MeshProcessor::ExtractBoundary3D(surface, mesh);
  return surface;
}

Mesh ExtractOpenSurface3D(const Mesh &boundary)
{
  Mesh surface;
  MeshProcessor::ExtractOpenSurface3D(surface, boundary);
  return surface;
}

Mesh MergeSurfaces3D(const std::vector<Mesh> &surfaces)
{
  Mesh merged_surface;
  MeshProcessor::MergeSurfaces3D(merged_surface, surfaces);
  return merged_surface;
}

} // namespace DTCC_BUILDER

PYBIND11_MODULE(_dtcc_builder, m)
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

  py::class_<DTCC_BUILDER::GridField>(m, "GridField")
      .def(py::init<>())
      .def_readonly("Grid", &DTCC_BUILDER::GridField::grid)
      .def_readonly("values", &DTCC_BUILDER::GridField::Values);

  py::class_<DTCC_BUILDER::Grid>(m, "Grid")
      .def(py::init<>())
      .def_readonly("xsize", &DTCC_BUILDER::Grid::XSize)
      .def_readonly("ysize", &DTCC_BUILDER::Grid::YSize)
      .def_readonly("xstep", &DTCC_BUILDER::Grid::XStep)
      .def_readonly("ystep", &DTCC_BUILDER::Grid::YStep);

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

  py::class_<DTCC_BUILDER::Mesh>(m, "Mesh")
      .def(py::init<>())
      .def_readonly("Vertices", &DTCC_BUILDER::Mesh::Vertices)
      .def_readonly("Faces", &DTCC_BUILDER::Mesh::Faces)
      .def_readonly("Normals", &DTCC_BUILDER::Mesh::Normals);

  py::class_<DTCC_BUILDER::VolumeMesh>(m, "VolumeMesh")
      .def(py::init<>())
      .def_readonly("numLayers", &DTCC_BUILDER::VolumeMesh::NumLayers)
      .def_readonly("Vertices", &DTCC_BUILDER::VolumeMesh::Vertices)
      .def_readonly("Cells", &DTCC_BUILDER::VolumeMesh::Cells)
      .def_readonly("Markers", &DTCC_BUILDER::VolumeMesh::Markers);

  m.def("createBuilderCityModel", &DTCC_BUILDER::createBuilderCityModel,
        "create builder point cloud from citymodel data");

  m.def("createBuilderPointCloud", &DTCC_BUILDER::createBuilderPointCloud,
        "create builder point cloud from numpy arrays");

  m.def("createBuilderGridField", &DTCC_BUILDER::createBuilderGridField,
        "create builder grid field from numpy arrays");

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
  m.def("GenerateVolumeMesh", &DTCC_BUILDER::GenerateVolumeMesh,
        "Generate 2D mesh");

  m.def("smooth_volume_mesh", &DTCC_BUILDER::smooth_volume_mesh,
        "Smooth volume mesh");
  m.def("TrimVolumeMesh", &DTCC_BUILDER::TrimVolumeMesh, "Trim 3D mesh");

  m.def("ExtractBoundary3D", &DTCC_BUILDER::ExtractBoundary3D,
        "Extract 3D boundary");

  m.def("GenerateSurface3D", &DTCC_BUILDER::GenerateSurfaces3D,
        "Generate 3D surface");

  m.def("MergeSurfaces3D", &DTCC_BUILDER::MergeSurfaces3D, "Merge 3D surfaces");
}