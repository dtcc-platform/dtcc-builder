#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

// Needs to come before JSON (nlohmann) include because of sloppy
// namespacing in VTK (typedef detail)...

//  DTCC includes

// #include "JSON.h"
#include <cmath>
#include <dolfin.h>
#include <iostream>

#include "LaplacianSmoother.h"

// Sandbox Includes
#include "include/JSON.hpp"
// #include "include/assembled.hpp"
// #include "include/fem.hpp"
#include "include/assembled.hpp"
#include "include/boundaryConditions.hpp"
#include "include/sparse.hpp"
#include "include/stiffnessMatrix.hpp"

using namespace DTCC;

void Help() { error("Usage: dtcc-sandbox Data_dir NumOfIterations"); }

void checkMeshFenics(Mesh3D &mesh,
                     CityModel &cityModel,
                     const GridField2D &dem,
                     double topHeight,
                     bool fixBuildings,
                     bool writeMat)
{
  info("Writing Stiffness Matrix A");

  // Disable dof reordering
  dolfin::parameters["reorder_dofs_serial"] = false;

  LaplacianSmoother::SmoothMesh3D(mesh, cityModel, dem, topHeight, fixBuildings,
                                  true);
  return;
}

void checkMeshGS(Mesh3D &mesh,
                 CityModel &cityModel,
                 const GridField2D &dem,
                 double topHeight,
                 bool fixBuildings,
                 bool writeMat)
{
  info("Writing Stiffness Matrix A: GS");

  const std::size_t nV = mesh.Vertices.size();

  stiffnessMatrix AK(mesh);

  std::vector<double> b(mesh.Vertices.size(), 0);

  BoundaryConditions bc(mesh, cityModel, dem, topHeight, fixBuildings);
  bc.apply(b);
  bc.apply(AK);

  std::cout << "Assembled A ... " << std::endl;
  std::vector<double> A = AK.assemble();

  for (size_t i = 0; i < nV; i++)
  {
    if (bc.vMarkers[i] > -4)
      A[nV * i + i] = 1;
  }

  std::stringstream ss;
  for (size_t i = 0; i < nV; i++)
  {
    for (size_t j = 0; j < nV; j++)
    {
      if (A[i * nV + j] != 0)
        ss << (i + 1) << " " << (j + 1) << " " << A[i * nV + j] << std::endl;
    }
  }

  std::ofstream fs;
  auto ff = "Matrix3_3.dat";
  if (fixBuildings)
  {
    ff = "Matrix3_5.dat";
  }
  fs.open(ff);
  fs << ss.rdbuf();
  fs.close();

  return;
}

int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (argc != 3)
  {
    Help();
    return 1;
  }
  std::string data_dir = argv[1];
  size_t max_iterations = atoi(argv[2]);

  Mesh3D mesh;
  CityModel cityModel;
  GridField2D dem;

  std::string mesh_dir = data_dir + "/CityMesh.json";
  std::string citymodel_dir = data_dir + "/CityModel.json";
  std::string dtm_dir = data_dir + "/DTM.json";

  JSON::Read(mesh, mesh_dir);
  JSON::Read(cityModel, citymodel_dir);
  JSON::Read(dem, dtm_dir);

  info(mesh);
  info(cityModel);
  info(dem);

  double topHeight{};
  topHeight = dem.Mean() + static_cast<double>(100.0);

  checkMeshFenics(mesh, cityModel, dem, topHeight, true, true);

  checkMeshGS(mesh, cityModel, dem, topHeight, true, true);

  return 0;
}
