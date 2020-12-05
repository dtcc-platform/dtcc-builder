// dtcc-sandbox
// Just a sandbox code for users to experiment with different API calls for
// core Vasilis Naserentin 2019
// Licensed under the MIT License

#include "CSV.h"
#include "CityJSON.h"
#include "JSON.h"
#include "Logging.h"
#include "VTK.h"

#include <iostream>
using namespace std;
using namespace DTCC;

void Help() { Error("Usage: dtcc-sandbox"); }

int main(int argc, char *argv[])
{
  if (argc > 4)
  {
    Help();
    return 1;
  }
  // std::string varname = argv[3];
  // std::string fileout = argv[2];
  // std::string filein = argv[1];
  // nlohmann::json json;
  // in.read_header(io::ignore_extra_column, "vendor", "size", "speed");
  // std::string vendor; int size; double speed;
  Mesh2D mesh2D;
  Mesh3D mesh3D;

  // JSON::Read(mesh2D, dataDirectory + "Mesh2D.json");
  JSON::Read(mesh2D, "/home/dtcc/core/build/Mesh2D.json");
  JSON::Read(mesh3D, "/home/dtcc/core/build/Mesh3D.json");

  std::cout << mesh2D.Cells[0] << std::endl;
  std::cout << mesh2D.Cells[1] << std::endl;
  std::cout << mesh2D.Cells.size() << std::endl;
  std::cout << mesh2D.Vertices.size() << std::endl;
  // std::cout<<mesh2D.Cells[0].size()<<std::endl;
  // std::cout<<mesh2D.Vertices[0].size()<<std::endl;

  /*
  Mesh3D mesh;
  mesh.Vertices.push_back(Point3D(0, 0, 0));
  mesh.Vertices.push_back(Point3D(1, 0, 0));
  mesh.Vertices.push_back(Point3D(0, 1, 0));
  mesh.Vertices.push_back(Point3D(0.5, 0.5, 1));
  mesh.Cells.push_back(Simplex3D(0, 1, 2, 3));
    std::cout<<mesh.Vertices[0].x<<std::endl;
  */
  // CSV csv;
  // csv.Read("test.csv", true);
  // CityJSON cityobj;
  // JSON::Read(cityobj,"test.json");
  // std::cout<<cityobj.CityObjects[0]<<std::endl;

  VTK::Write2(mesh2D, "TEST2D.vtu");
  VTK::Write3(mesh3D, "TEST3D.vtu");

  //DTCC::CityJSON::Read("HI");
}
