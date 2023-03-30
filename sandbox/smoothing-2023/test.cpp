#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <chrono>

// Needs to come before JSON (nlohmann) include because of sloppy
// namespacing in VTK (typedef detail)...

// #include "fem.h"
//  DTCC includes
//#include "JSON.h"
#include "Mesh.h"
#include "Timer.h"

//#include "./include/fem.hpp"
//#include "./include/sparse.hpp"
//#include "./include/assembled.hpp"

/// Read JSON data from file
/* static void Read(nlohmann::json& json, std::string fileName)
{
  std::cout<<"JSON: Reading from file " + fileName + "..."<<std::endl;
  std::ifstream f(fileName);
  if (!f)
    error("Unable to read from file " + fileName);
  f >> json;
} */


using namespace DTCC_BUILDER;

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

void Help() { error("Usage: dtcc-sandbox CityMesh.json numOfIterations"); }
void test_stiffnessMatrix(Mesh3D * mesh);

int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (argc != 3)
  {
    Help();
    return 1;
  }

  // Read parameters
  //Parameters p;
  Mesh3D m;
  //Read(m, argv[1]);

  size_t max_iterations = atoi(argv[2]);
  std::cout<< "Max iter "<< max_iterations <<std::endl;
  info(m);
  
  std::cout << "\nNumber of Cells: " << m.Cells.size() << std::endl;
  std::cout << "Number of Vertices: " << m.Vertices.size() << std::endl;
  return 0;
  }