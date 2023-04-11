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
#include "include/JSON.hpp"
// #include "JSON.h"
#include "LaplacianSmoother.h"
#include "Mesh.h"
#include "Timer.h"
#include "datamodel/CityModel.h"

#include "include/assembled.hpp"
#include "include/fem.hpp"
#include "include/sparse.hpp"

using namespace DTCC;

void Help() { error("Usage: dtcc-sandbox Data_dir NumOfIterations"); }

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

  Mesh3D m;
  CityModel cm;
  GridField2D dtm;

  std::string mesh_dir = data_dir + "/CityMesh.json";
  std::string citymodel_dir = data_dir + "/CityModel.json";
  std::string dtm_dir = data_dir + "/DTM.json";

  JSON::Read(m, mesh_dir);
  JSON::Read(cm, citymodel_dir);
  JSON::Read(dtm, dtm_dir);

  std::cout << "Max iterations for Iterative Solver: " << max_iterations
            << std::endl;
  info(m);
  info(cm);
  info(dtm);
  std::cout << "\nNumber of Cells: " << m.Cells.size() << std::endl;
  std::cout << "Number of Vertices: " << m.Vertices.size() << std::endl;
  std::cout << "Number of Markers: " << m.Markers.size() << std::endl;
  std::cout << "Number of Buildings: " << cm.Buildings.size() << std::endl;

  int *vMarkers = new int[m.Vertices.size()];
  getVerticeMarkers(m, vMarkers);

  double topHeight{};
  topHeight = dtm.Mean() + static_cast<double>(100.0);
  {
    Timer timer("Fenics Smooth 3D mesh");
    topHeight = dtm.Mean() + static_cast<double>(100.0);
    LaplacianSmoother::SmoothMesh3D(m, cm, dtm, topHeight, false, false);
    timer.Stop();
    timer.Print();
  }

  // {
  //   Timer timer("Jacobi Smooth 3D mesh");
  //   smoothLaplaceJacobi(m, cm, dtm, max_iterations);
  //   timer.Stop();
  //   timer.Print();
  // }

  // {
  //   Timer timer("Jacobi Smooth 3D mesh");
  //   // assembled_GaussSeidel(&m, &cm, &dtm, max_iterations);
  //   timer.Stop();
  //   timer.Print();
  // }

  delete[] vMarkers;
  return 0;
}