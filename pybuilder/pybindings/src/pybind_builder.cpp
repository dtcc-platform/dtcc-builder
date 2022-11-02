#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "datamodel/CityModel.h"
// DTCC includes
#include "CityModelGenerator.h"
#include "CommandLine.h"
#include "ElevationModelGenerator.h"
#include "GridField.h"
#include "JSON.h"
#include "LAS.h"
#include "Logging.h"
#include "ParameterProcessor.h"
#include "Parameters.h"
#include "PointCloud.h"
#include "Polygon.h"
#include "SHP.h"
#include "Timer.h"
#include "Utils.h"
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

PointCloud GlobalOutlierRemover(PointCloud &pointCloud, double outlierMargin)
{
  PointCloudProcessor::RemoveOutliers(pointCloud, outlierMargin);
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

} // namespace DTCC_BUILDER

int add(int i, int j) { return i + j; }

PYBIND11_MODULE(_pybuilder, m)
{

  py::class_<DTCC_BUILDER::CityModel>(m, "CityModel").def(py::init<>());

  py::class_<DTCC_BUILDER::PointCloud>(m, "PointCloud").def(py::init<>());

  py::class_<DTCC_BUILDER::GridField2D>(m, "GridField2D").def(py::init<>());

  m.doc() = "python bindings for dtcc-builder";

  m.def("add", &add, "A function that adds two numbers");

  m.def("GenerateCityModel", &DTCC_BUILDER::GenerateCityModel,
        "load shp file into city model");

  m.def("LASReadDirectory", &DTCC_BUILDER::LASReadDirectory,
        "load all .las files in directory");

  m.def("LASReadFile", &DTCC_BUILDER::LASReadFile, "load .las file");

  m.def("LASBounds", &DTCC_BUILDER::LASBounds,
        "calculate bounding box of all .las files in directorty");

  m.def("GenerateElevationModel", &DTCC_BUILDER::GenerateElevationModel,
        "generate height field from point cloud");

  m.def("GlobalOutlierRemover", &DTCC_BUILDER::GlobalOutlierRemover,
        "Remove all points more than a given number of standard deviations "
        "from the mean for the z-coordinate");
}
