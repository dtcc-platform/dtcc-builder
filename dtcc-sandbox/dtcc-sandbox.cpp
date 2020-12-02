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
  Mesh3D mesh;

  CSV csv;
  //csv.Read("test.csv", true);
  CityJSON cityobj;
  // JSON::Read(cityobj,"test.json");
  // std::cout<<cityobj.CityObjects[0]<<std::endl;

  std::string filename = "test";

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->InsertNextPoint(0, 0, 0);
  points->InsertNextPoint(1, 0, 0);
  points->InsertNextPoint(1, 1, 0);
  //  points->InsertNextPoint(0, 1, 1);

  vtkSmartPointer<vtkTriangle> tetra = vtkSmartPointer<vtkTriangle>::New();

  tetra->GetPointIds()->SetId(0, 0);
  tetra->GetPointIds()->SetId(1, 1);
  tetra->GetPointIds()->SetId(2, 2);
  //  tetra->GetPointIds()->SetId(3, 3);

  vtkSmartPointer<vtkCellArray> cellArray =
      vtkSmartPointer<vtkCellArray>::New();
  cellArray->InsertNextCell(tetra);

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
