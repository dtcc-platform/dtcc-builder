#include "JSON.h"
#include "Mesh.h"
#include "MeshField.h"
#include "VTK.h"
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>

using namespace DTCCBUILDER;

int main(int, char *[])
{
  GridField2D dtm; //, dtm;
  JSON::Read(dtm, "/home/dtcc/core/data/DTM.json");
  // JSON::Read(dtm, "/home/dtcc/core/data/DTM.json");
  std::cout << dtm.Grid.XSize << std::endl;
  std::cout << dtm.Grid.YSize << std::endl;
  std::cout << dtm.Grid.XStep << std::endl;
  std::cout << dtm.Grid.YStep << std::endl;
  std::cout << dtm.Values.size() << std::endl;

  vtkSmartPointer<vtkStructuredGrid> structuredGrid =
      vtkSmartPointer<vtkStructuredGrid>::New();

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  // structuredGrid->SetDimensions(dtm.Grid.XSize, dtm.Grid.YSize,0);
  structuredGrid->SetDimensions(dtm.Grid.XSize, dtm.Grid.YSize, 1);

  for (uint i = 0; i < dtm.Grid.XSize; ++i)
    for (uint j = 0; j < dtm.Grid.YSize; ++j)
    {
      {
        points->InsertNextPoint(i * dtm.Grid.XStep, j * dtm.Grid.YStep, 0);
      }
    }

  std::cout << "Done doing the points" << std::endl;

  structuredGrid->SetPoints(points);
  std::cout << "Done copying the points" << std::endl;

  vtkFloatArray *Values = vtkFloatArray::New();
  Values->SetName("Values");
  for (uint i = 0; i < dtm.Values.size(); i++)
    Values->InsertTuple1(i, dtm.Values[i]);
  structuredGrid->GetPointData()->SetScalars(Values);
  std::cout << "Done with tuples" << std::endl;

  // Write file
  vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
      vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
  writer->SetFileName("output2.vts");
  writer->SetInputData(structuredGrid);
  writer->Write();

  VTK::Write(dtm, "output_VTK.vts");

  return EXIT_SUCCESS;
}
