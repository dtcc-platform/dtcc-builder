// dtcc-sandbox
// Just a sandbox code for users to experiment with different API calls for
// core Vasilis Naserentin 2019
// Licensed under the MIT License

#include "CSV.h"
#include "CityJSON.h"
#include "JSON.h"
#include "Logging.h"
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkXMLUnstructuredGridReader.h>
//#include <vtkDataSetMapper.h>
//#include <vtkActor.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

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

  std::string filename = "test2D.vtu";

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  for (uint i = 0; i < mesh2D.Vertices.size(); i++)
  {
    points->InsertNextPoint(mesh2D.Vertices[i].x, mesh2D.Vertices[i].y, 0);
  }

  vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

  vtkSmartPointer<vtkCellArray> cellArray =
      vtkSmartPointer<vtkCellArray>::New();

  // vtkSmartPointer<vtkTriangle> triangle2 =
  // vtkSmartPointer<vtkTriangle>::New();
  for (uint i = 0; i < mesh2D.Cells.size(); i++)
  // for (uint i = 0; i < 100; i++)

  {
    triangle->GetPointIds()->SetId(0, mesh2D.Cells[i].v0);
    triangle->GetPointIds()->SetId(1, mesh2D.Cells[i].v1);
    triangle->GetPointIds()->SetId(2, mesh2D.Cells[i].v2);

    cellArray->InsertNextCell(triangle);
    //  triangle->Reset();
  }

  //  triangle2->GetPointIds()->SetId(0, mesh2D.Cells[1].v0);
  //  triangle2->GetPointIds()->SetId(1, mesh2D.Cells[1].v1);
  //  triangle2->GetPointIds()->SetId(2, mesh2D.Cells[1].v2);

  /*
    for (uint i = 0; i<3;i++)
    {
    points->InsertNextPoint(mesh2D.Vertices[i].x, mesh2D.Vertices[i].y,0);
    }
  */
  // Find first triangle
  //   points->InsertNextPoint(mesh2D.Vertices[0].x, mesh2D.Vertices[0].y, 0);
  ///   points->InsertNextPoint(mesh2D.Vertices[1].x, mesh2D.Vertices[1].y, 0);
  //  points->InsertNextPoint(mesh2D.Vertices[2].x, mesh2D.Vertices[2].y, 0);

  // points->InsertNextPoint(1, 0, 0);
  // points->InsertNextPoint(1, 1, 0);
  // points->InsertNextPoint(0, 1, 1);


  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
  unstructuredGrid->SetPoints(points);
  unstructuredGrid->SetCells(VTK_TRIANGLE, cellArray);

  // Write file
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
      vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData(unstructuredGrid);
  writer->Write();

  //DTCC::CityJSON::Read("HI");
}
