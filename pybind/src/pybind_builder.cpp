#include <pybind11/pybind11.h>

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
#include "Polygon.h"
#include "SHP.h"
#include "Timer.h"
#include "Utils.h"
#include "VertexSmoother.h"

namespace py = pybind11;

namespace DTCC_BUILDER
{

CityModel GenerateCityModel(std::string shp_file, double minBuildingDistance,double minBuildingSize)
{
  std::vector<Polygon> footprints;
  std::vector<std::string> UUIDs;
  std::vector<int> entityIDs;

  auto loading_timer = Timer("load data");
  SHP::Read(footprints, shp_file, &UUIDs, &entityIDs);
  info("Loaded " + str(footprints.size()) + " building footprints");

  CityModel cityModel;

  return cityModel;



}
} // namespace DTCC_BUILDER


int add(int i, int j) {
    return i + j;
}


PYBIND11_MODULE(_pybuilder, m) {

    py::class_<DTCC_BUILDER::CityModel>(m,"CityModel")
        .def(py::init<>());

    m.doc() = "python bindings for dtcc-builder"; 

    m.def("add", &add, "A function that adds two numbers");

    m.def("GenerateCityModel", &DTCC_BUILDER::GenerateCityModel, 
      "load shp file into city model");
}
