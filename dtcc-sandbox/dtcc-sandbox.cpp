#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>

int main(int, char *[])
{
  vtkSmartPointer<vtkStructuredGrid> structuredGrid =
      vtkSmartPointer<vtkStructuredGrid>::New();

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
      {
        {
          {
            points->InsertNextPoint(i, j, k);
          }
        }
      }

  //
  structuredGrid->SetDimensions(3, 3, 3);
  structuredGrid->SetPoints(points);

  vtkFloatArray *scalars = vtkFloatArray::New();
  for (int i = 0; i < 27; i++)
    scalars->InsertTuple1(i, 1.0 * i / (27));
  structuredGrid->GetPointData()->SetScalars(scalars);

  // Write file
  vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
      vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
  writer->SetFileName("output.vts");
  writer->SetInputData(structuredGrid);
  writer->Write();

  return EXIT_SUCCESS;
}
