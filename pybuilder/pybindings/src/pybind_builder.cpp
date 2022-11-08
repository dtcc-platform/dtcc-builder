#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "datamodel/Building.h"
#include "datamodel/CityModel.h"

// Needs to come before JSON (nlohmann) include because of sloppy
// namespacing in VTK (typedef detail)...
#include "VTK.h"

// DTCC includes
#include "CityModelGenerator.h"
#include "ElevationModelGenerator.h"
#include "GridField.h"
#include "JSON.h"
#include "LAS.h"
#include "LaplacianSmoother.h"
#include "Logging.h"
#include "Mesh.h"
#include "MeshGenerator.h"
#include "MeshIO.h"
#include "MeshProcessor.h"
#include "Point.h"
#include "PointCloud.h"
#include "Polygon.h"
#include "SHP.h"
#include "Surface.h"
#include "VertexSmoother.h"

namespace py = pybind11;

namespace DTCC_BUILDER
{

// CityModel
CityModel GenerateCityModel(std::string shp_file,
                            py::tuple bounds,
                            double minBuildingDistance,
                            double minBuildingSize)
{
  std::vector<Polygon> footprints;
  std::vector<std::string> UUIDs;
  std::vector<int> entityIDs;

  auto loading_timer = Timer("load data");
  SHP::Read(footprints, shp_file, &UUIDs, &entityIDs);
  info("Loaded " + str(footprints.size()) + " building footprints");

  CityModel cityModel;
  double px = bounds[0].cast<double>();
  double py = bounds[1].cast<double>();
  double qx = bounds[2].cast<double>();
  double qy = bounds[3].cast<double>();
  auto bbox = BoundingBox2D(Point2D(px, py), Point2D(qx, qy));

  CityModelGenerator::GenerateCityModel(cityModel, footprints, UUIDs, entityIDs,
                                        bbox, minBuildingDistance,
                                        minBuildingSize);

  return cityModel;
}

CityModel SetCityModelOrigin(CityModel &cityModel, py::tuple origin)
{
  double px = origin[0].cast<double>();
  double py = origin[1].cast<double>();
  Point2D o = Point2D(px, py);
  cityModel.SetOrigin(o);

  return cityModel;
}

CityModel CleanCityModel(CityModel &cityModel, double minVertDistance)
{
  CityModelGenerator::CleanCityModel(cityModel, minVertDistance);
  return cityModel;
}

CityModel ExtractBuildingPoints(CityModel &cityModel,
                                const PointCloud &pointCloud,
                                double groundMargin,
                                double groundOutlierMargin)
{
  CityModelGenerator::ExtractBuildingPoints(cityModel, pointCloud, groundMargin,
                                            groundOutlierMargin);
  return cityModel;
}

CityModel BuildingPointsRANSACOutlierRemover(CityModel &cityModel,
                                             double outlier_margin,
                                             size_t iterations)
{
  CityModelGenerator::BuildingPointsRANSACOutlierRemover(
      cityModel, outlier_margin, iterations);
  return cityModel;
}

CityModel BuildingPointsOutlierRemover(CityModel &cityModel,
                                       size_t neighbors,
                                       double outlier_margin)
{
  CityModelGenerator::BuildingPointsOutlierRemover(cityModel, neighbors,
                                                   outlier_margin);
  return cityModel;
}

CityModel ComputeBuildingHeights(CityModel &cityModel,
                                 const GridField2D &dtm,
                                 double groundPercentile,
                                 double roofPercentile)
{
  CityModelGenerator::ComputeBuildingHeights(cityModel, dtm, groundPercentile,
                                             roofPercentile);
  return cityModel;
}

void WriteCityModelJSON(const CityModel &cityModel, std::string path)
{
  JSON::Write(cityModel, path, cityModel.Origin);
}

CityModel ReadCityModelJSON(std::string path)
{
  CityModel cityModel;
  JSON::Read(cityModel, path);
  return cityModel;
}

// PointCloud
PointCloud LASReadDirectory(std::string las_directory, bool extra_data = true)
{
  PointCloud pc;
  LAS::ReadDirectory(pc, las_directory, extra_data);
  return pc;
}

PointCloud LASReadFile(std::string las_file, bool extra_data = true)
{
  PointCloud pc;
  LAS::Read(pc, las_file, extra_data);
  return pc;
}

py::tuple LASBounds(std::string las_directory)
{
  BoundingBox2D bb;
  LAS::BoundsDirectory(bb, las_directory);
  py::tuple bbox = py::make_tuple(bb.P.x, bb.P.y, bb.Q.x, bb.Q.y);
  return bbox;
}

PointCloud SetPointCloudOrigin(PointCloud &pointCloud, py::tuple origin)
{
  double px = origin[0].cast<double>();
  double py = origin[1].cast<double>();
  Point2D o = Point2D(px, py);
  pointCloud.SetOrigin(o);
  return pointCloud;
}

PointCloud GlobalOutlierRemover(PointCloud &pointCloud, double outlierMargin)
{
  PointCloudProcessor::RemoveOutliers(pointCloud, outlierMargin);
  return pointCloud;
}

PointCloud VegetationFilter(PointCloud &pointCloud)
{
  PointCloudProcessor::NaiveVegetationFilter(pointCloud);
  return pointCloud;
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
                    double topHeight,
                    bool fixBuildings)
{
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

bool WriteSurface3D(const Surface3D &surface,
                    std::string fileName,
                    std::string format,
                    bool YUp)
{
  MeshIO::Write(surface, fileName, format, YUp);
  return true;
}

bool WriteVTKMesh3D(const Mesh3D &mesh, std::string filepath)
{
  VTK::Write(mesh, filepath);
  return true;
}

bool WriteVTKMesh2D(const Mesh2D &mesh, std::string filepath)
{
  VTK::Write(mesh, filepath);
  return true;
}

} // namespace DTCC_BUILDER

PYBIND11_MODULE(_pybuilder, m)
{

  py::class_<DTCC_BUILDER::CityModel>(m, "CityModel")
      .def(py::init<>())
      .def(
          "__len__",
          [](const DTCC_BUILDER::CityModel &cm) { return cm.Buildings.size(); })
      .def_readonly("buildings", &DTCC_BUILDER::CityModel::Buildings);

  py::class_<DTCC_BUILDER::Building>(m, "Building")
      .def(py::init<>())
      .def_readwrite("error", &DTCC_BUILDER::Building::error)
      .def_readwrite("uuid", &DTCC_BUILDER::Building::UUID)
      .def_readwrite("propertyUUID", &DTCC_BUILDER::Building::PropertyUUID)
      .def_readwrite("height", &DTCC_BUILDER::Building::Height)
      .def_readwrite("groundHeight", &DTCC_BUILDER::Building::GroundHeight)
      .def_readonly("footprint", &DTCC_BUILDER::Building::Footprint)
      .def_readonly("grounPoints", &DTCC_BUILDER::Building::GroundPoints)
      .def_readonly("roofPoints", &DTCC_BUILDER::Building::RoofPoints);

  py::class_<DTCC_BUILDER::Point2D>(m, "Point2D")
      .def(py::init<>())
      .def("__repr__",
           [](const DTCC_BUILDER::Point3D &p) {
             return "<Point3D (" + DTCC_BUILDER::str(p.x) + ", " +
                    DTCC_BUILDER::str(p.y) + ")>";
           })
      .def_readonly("x", &DTCC_BUILDER::Point2D::x)
      .def_readonly("y", &DTCC_BUILDER::Point2D::y);

  py::class_<DTCC_BUILDER::Point3D>(m, "Point3D")
      .def(py::init<>())
      .def("__repr__",
           [](const DTCC_BUILDER::Point3D &p) {
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
                    &DTCC_BUILDER::PointCloud::Classifications);

  py::class_<DTCC_BUILDER::Simplex2D>(m, "Simplex2D").def(py::init<>());

  py::class_<DTCC_BUILDER::Simplex3D>(m, "Simplex3D").def(py::init<>());

  py::class_<DTCC_BUILDER::Grid2D>(m, "Grid2D")
      .def(py::init<>())
      .def_readonly("BoundingBox", &DTCC_BUILDER::Grid2D::BoundingBox)
      .def_readonly("XSize", &DTCC_BUILDER::Grid2D::XSize)
      .def_readonly("YSize", &DTCC_BUILDER::Grid2D::YSize)
      .def_readonly("XStep", &DTCC_BUILDER::Grid2D::XStep)
      .def_readonly("XStep", &DTCC_BUILDER::Grid2D::YStep);

  py::class_<DTCC_BUILDER::GridField2D>(m, "GridField2D")
      .def(py::init<>())
      .def_readonly("Grid", &DTCC_BUILDER::GridField2D::Grid);

  py::class_<DTCC_BUILDER::Mesh2D>(m, "Mesh2D").def(py::init<>());

  py::class_<DTCC_BUILDER::Mesh3D>(m, "Mesh3D")
      .def(py::init<>())
      .def_readonly("numLayers", &DTCC_BUILDER::Mesh3D::NumLayers);

  py::class_<DTCC_BUILDER::Surface3D>(m, "Surface3D").def(py::init<>());

  m.doc() = "python bindings for dtcc-builder";

  m.def("GenerateCityModel", &DTCC_BUILDER::GenerateCityModel,
        "load shp file into city model");

  m.def("SetCityModelOrigin", &DTCC_BUILDER::SetCityModelOrigin,
        "Set Origin on CityModel");

  m.def("CleanCityModel", &DTCC_BUILDER::CleanCityModel,
        "clean city model polygons");

  m.def("ExtractBuildingPoints", &DTCC_BUILDER::ExtractBuildingPoints,
        "extarct points from point cloud for each building");

  m.def("LASReadDirectory", &DTCC_BUILDER::LASReadDirectory,
        "load all .las files in directory");

  m.def("LASReadFile", &DTCC_BUILDER::LASReadFile, "load .las file");

  m.def("LASBounds", &DTCC_BUILDER::LASBounds,
        "calculate bounding box of all .las files in directorty");

  m.def("SetPointCloudOrigin", &DTCC_BUILDER::SetPointCloudOrigin,
        "set point cloud origin");

  m.def("GlobalOutlierRemover", &DTCC_BUILDER::GlobalOutlierRemover,
        "Remove all points more than a given number of standard deviations "
        "from the mean for the z-coordinate");

  m.def("VegetationFilter", &DTCC_BUILDER::VegetationFilter,
        "Remove possible vegetation filters");

  m.def(
      "BuildingPointsRANSACOutlierRemover",
      &DTCC_BUILDER::BuildingPointsRANSACOutlierRemover,
      "Use RANSAC to remove extreme outliers. Only useful on very noisy data");

  m.def("BuildingPointsOutlierRemover",
        &DTCC_BUILDER::BuildingPointsOutlierRemover,
        "remove outliers from roof points using statistcal outlier algorithm");

  m.def("ComputeBuildingHeights", &DTCC_BUILDER::ComputeBuildingHeights,
        "Calculate building heights based on point cloud and dtm data");

  m.def("WriteCityModelJSON", &DTCC_BUILDER::WriteCityModelJSON,
        "Write CityModel to JSON");

  m.def("ReadCityModelJSON", &DTCC_BUILDER::ReadCityModelJSON,
        "Load CityModel from JSON");

  m.def("GenerateElevationModel", &DTCC_BUILDER::GenerateElevationModel,
        "generate height field from point cloud");

  m.def("SmoothElevation", &DTCC_BUILDER::SmoothElevation,
        "Smooth  elevation grid field");

  m.def("GenerateMesh2D", &DTCC_BUILDER::GenerateMesh2D, "Generate 2D mesh");

  m.def("GenerateMesh3D", &DTCC_BUILDER::GenerateMesh3D, "Generate 3D mesh");

  m.def("SmoothMesh3D", &DTCC_BUILDER::SmoothMesh3D,
        "Laplacian smooth 3D mesh");

  m.def("TrimMesh3D", &DTCC_BUILDER::TrimMesh3D,
        "Trim 3D mesh. The mesh is trimmed by removing tetrahedra inside "
        "building shape.");

  m.def("GenerateSurfaces3D", &DTCC_BUILDER::GenerateSurfaces3D,
        "Generate 3D surface meshes for visualization. Returns a list of "
        "surfaces. The first surface is the ground followed buy each building");

  m.def("ExtractBoundary3D", &DTCC_BUILDER::ExtractBoundary3D,
        "Extract the boundary of a 3D mesh as a 3D surface.");

  m.def("ExtractOpenSurface3D", &DTCC_BUILDER::ExtractOpenSurface3D,
        "Extract an open surface from a boundary, excluding top and sides");

  m.def("MergeSurfaces3D", &DTCC_BUILDER::MergeSurfaces3D,
        "Merge surfaces into a single surface");

  m.def("WriteSurface3D", &DTCC_BUILDER::WriteSurface3D,
        "Write surface to file");

  m.def("WriteVTKMesh3D", &DTCC_BUILDER::WriteVTKMesh3D,
        "Write 3D mesh to VTK format");

  m.def("WriteVTKMesh2D", &DTCC_BUILDER::WriteVTKMesh2D,
        "Write 3D mesh to VTK format");
}
