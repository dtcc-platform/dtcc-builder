// Copyright (C) 2023 Dag Wästberg
// Licensed under the MIT License
//
// Modified by Anders Logg 2023

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <pybind11/stl.h>

#include "CityBuilder.h"
#include "ElevationBuilder.h"
#include "MeshBuilder.h"
#include "MeshProcessor.h"
#include "Smoother.h"
#include "VertexSmoother.h"
#include "model/Building.h"
#include "model/City.h"
#include "model/GridField.h"
#include "model/Mesh.h"
#include "model/Point.h"
#include "model/PointCloud.h"
#include "model/Polygon.h"

namespace py = pybind11;

namespace DTCC_BUILDER
{

Polygon create_polygon(py::list vertices, py::list holes)
{
  Polygon poly;
  for (size_t i = 0; i < vertices.size(); i++)
  {
    auto pt = vertices[i].cast<py::tuple>();
    poly.Vertices.push_back(
        Point2D(pt[0].cast<double>(), pt[1].cast<double>()));
  }
  if (holes.size() > 0)
  {
    for (size_t i = 0; i < holes.size(); i++)
    {
      auto hl = holes[i].cast<py::list>();
      std::vector<Point2D> hole;
      for (size_t j = 0; j < hl.size(); j++)
      {
        auto pt = hl[j].cast<py::tuple>();
        hole.push_back(Point2D(pt[0].cast<double>(), pt[1].cast<double>()));
      }
      poly.Holes.push_back(hole);
    }
  }
  return poly;
}

City create_city(py::list footprints,
                 py::list holes,
                 py::list uuids,
                 py::list heights,
                 py::list ground_levels,
                 py::tuple origin)
{
  City city;

  city.Origin = Point2D(origin[0].cast<double>(), origin[1].cast<double>());
  size_t num_buildings = footprints.size();
  for (size_t i = 0; i < num_buildings; i++)
  {
    auto footprint = footprints[i].cast<py::list>();
    auto uuid = uuids[i].cast<std::string>();
    auto height = heights[i].cast<double>();
    auto ground_level = ground_levels[i].cast<double>();
    auto hls = holes[i].cast<py::list>();

    Polygon poly;
    for (size_t j = 0; j < footprint.size(); j++)
    {
      auto pt = footprint[j].cast<py::tuple>();

      poly.Vertices.push_back(
          Point2D(pt[0].cast<double>(), pt[1].cast<double>()));
    }
    if (hls.size() > 0)
    {
      for (size_t j = 0; j < hls.size(); j++)
      {
        auto hl = hls[j].cast<py::list>();
        std::vector<Point2D> hole;
        for (size_t k = 0; k < hl.size(); k++)
        {
          auto pt = hl[k].cast<py::tuple>();
          hole.push_back(Point2D(pt[0].cast<double>(), pt[1].cast<double>()));
        }
        poly.Holes.push_back(hole);
      }
    }
    Building building;
    building.Footprint = poly;
    building.UUID = uuid;
    building.Height = height;
    building.GroundHeight = ground_level;
    city.Buildings.push_back(building);
  }
  // Cleaning should be done elsewhere (in Python)
  // This function should only convert
  // CityBuilder::clean_city(city, 1.0);

  return city;
}

PointCloud create_pointcloud(py::array_t<double> pts,
                             py::array_t<uint8_t> cls,
                             py::array_t<uint8_t> retNumber,
                             py::array_t<uint8_t> numReturns)
{
  auto pts_r = pts.unchecked<2>();
  auto cls_r = cls.unchecked<1>();
  auto retNumber_r = retNumber.unchecked<1>();
  auto numReturns_r = numReturns.unchecked<1>();

  auto pt_count = pts_r.shape(0);

  bool hasClassification = cls_r.size() == pt_count;
  bool hasReturnNumber = retNumber_r.size() == pt_count;
  bool hasNumberOfReturns = numReturns_r.size() == pt_count;

  PointCloud pointCloud;
  for (int i = 0; i < pt_count; i++)
  {
    pointCloud.Points.push_back(Point3D(pts_r(i, 0), pts_r(i, 1), pts_r(i, 2)));
    if (hasClassification)
      pointCloud.Classifications.push_back(cls_r(i));
    else
      pointCloud.Classifications.push_back(1);
    if (hasReturnNumber and hasNumberOfReturns)
      pointCloud.ScanFlags.push_back(
          PointCloudProcessor::packScanFlag(retNumber_r(i), numReturns_r(i)));
  }

  pointCloud.BuildHasClassifications();
  pointCloud.CalculateBoundingBox();
  return pointCloud;
}

GridField create_gridfield(py::array_t<double> data,
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

} // namespace DTCC_BUILDER

PYBIND11_MODULE(_dtcc_builder, m)
{
  py::class_<DTCC_BUILDER::City>(m, "City")
      .def(py::init<>())
      .def("__len__",
           [](const DTCC_BUILDER::City &cm) { return cm.Buildings.size(); })
      .def_readonly("buildings", &DTCC_BUILDER::City::Buildings)
      .def_readonly("origin", &DTCC_BUILDER::City::Origin);

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
      .def_readonly("vertices", &DTCC_BUILDER::Polygon::Vertices)
      .def_readonly("holes", &DTCC_BUILDER::Polygon::Holes);

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
      .def_readonly("num_layers", &DTCC_BUILDER::VolumeMesh::num_layers)
      .def_readonly("Vertices", &DTCC_BUILDER::VolumeMesh::Vertices)
      .def_readonly("Cells", &DTCC_BUILDER::VolumeMesh::Cells)
      .def_readonly("Markers", &DTCC_BUILDER::VolumeMesh::Markers);

  m.def("create_polygon", &DTCC_BUILDER::create_polygon, "Create C++ polygon");

  m.def("create_city", &DTCC_BUILDER::create_city, "Create C++ city");

  m.def("create_pointcloud", &DTCC_BUILDER::create_pointcloud,
        "Create C++ point cloud");

  m.def("create_gridfield", &DTCC_BUILDER::create_gridfield,
        "Create C++ grid field");

  m.def("remove_vegetation",
        &DTCC_BUILDER::PointCloudProcessor::remove_vegetation,
        "Remove vegetation from point cloud");

  m.def("remove_building_point_outliers_statistical",
        &DTCC_BUILDER::CityBuilder::remove_building_point_outliers_statistical,
        "Remove building point outliers (statistical)");

  m.def("remove_building_point_outliers_ransac",
        &DTCC_BUILDER::CityBuilder::remove_building_point_outliers_ransac,
        "Remove building point outliers (RANSAC)");

  m.def("compute_building_points",
        &DTCC_BUILDER::CityBuilder::compute_building_points,
        "Compute building points from point cloud");

  m.def("build_elevation", &DTCC_BUILDER::ElevationBuilder::build_elevation,
        "Build height field from point cloud");

  m.def("smooth_field", &DTCC_BUILDER::VertexSmoother::smooth_field,
        "Smooth grid field");

  m.def("clean_city", &DTCC_BUILDER::CityBuilder::clean_city, "Clean city");

  m.def("build_mesh", &DTCC_BUILDER::MeshBuilder::build_mesh,
        "Build mesh for city, returning a list of meshes");

  m.def("build_ground_mesh", &DTCC_BUILDER::MeshBuilder::build_ground_mesh,
        "Build ground mesh");

  m.def("layer_ground_mesh", &DTCC_BUILDER::MeshBuilder::layer_ground_mesh,
        "Layer ground mesh");

  m.def("smooth_volume_mesh", &DTCC_BUILDER::Smoother::smooth_volume_mesh,
        "Smooth volume mesh");

  m.def("trim_volume_mesh", &DTCC_BUILDER::MeshBuilder::trim_volume_mesh,
        "Trim volume mesh by removing cells inside buildings");

  m.def("extrude_footprint", &DTCC_BUILDER::MeshBuilder::extrude_footprint,
        "Extrude footprint to a mesh");

  m.def("compute_boundary_mesh",
        &DTCC_BUILDER::MeshProcessor::compute_boundary_mesh,
        "Compute boundary mesh from volume mesh");

  m.def("compute_open_mesh", &DTCC_BUILDER::MeshProcessor::compute_open_mesh,
        "Compute open mesh from boundary, excluding top and sides");

  m.def("merge_meshes", &DTCC_BUILDER::MeshProcessor::merge_meshes,
        "Merge meshes into a single mesh");
}
