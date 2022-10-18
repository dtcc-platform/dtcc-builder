// Copyright (C) 2020 Anders Logg and Vasilis Naserentin
// Licensed under the MIT License

#ifndef DTCC_VTK_H
#define DTCC_VTK_H

#include <fstream>
#include <iostream>
#include <vector>

#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "GridField.h"
#include "Logging.h"
#include "Mesh.h"
#include "Simplex.h"
#include "Surface.h"
#include "Vector.h"

namespace DTCC
{

class VTK
{
public:
  /// Write 3D mesh to VTK (.vtu) file.
  ///
  /// @param mesh3D The 3D mesh
  /// @parame fileName Filename (path)
  static void Write(const Mesh3D &mesh3D, std::string fileName)
  {
    info("VTK: Writing 3D mesh to file " + fileName);

    // Set vertex data
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (uint i = 0; i < mesh3D.Vertices.size(); i++)
    {
      points->InsertNextPoint(mesh3D.Vertices[i].x, mesh3D.Vertices[i].y,
                              mesh3D.Vertices[i].z);
    }

    // Set cell data
    vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();
    vtkSmartPointer<vtkCellArray> cellArray =
        vtkSmartPointer<vtkCellArray>::New();
    for (uint i = 0; i < mesh3D.Cells.size(); i++)

    {
      tetra->GetPointIds()->SetId(0, mesh3D.Cells[i].v0);
      tetra->GetPointIds()->SetId(1, mesh3D.Cells[i].v1);
      tetra->GetPointIds()->SetId(2, mesh3D.Cells[i].v2);
      tetra->GetPointIds()->SetId(3, mesh3D.Cells[i].v3);
      cellArray->InsertNextCell(tetra);
    }

    // Set grid data
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_TETRA, cellArray);

    // Write file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
  }

  /// Write 2D mesh to VTK (.vtu) file.
  ///
  /// @param mesh2D The 2D mesh
  /// @parame fileName Filename (path)
  static void Write(const Mesh2D &mesh2D, std::string fileName)
  {
    info("VTK: Writing 2D mesh to file " + fileName);

    // Set vertex data
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (uint i = 0; i < mesh2D.Vertices.size(); i++)
    {
      points->InsertNextPoint(mesh2D.Vertices[i].x, mesh2D.Vertices[i].y, 0);
    }

    // Set cell data
    vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
    vtkSmartPointer<vtkCellArray> cellArray =
        vtkSmartPointer<vtkCellArray>::New();
    for (uint i = 0; i < mesh2D.Cells.size(); i++)
    {
      triangle->GetPointIds()->SetId(0, mesh2D.Cells[i].v0);
      triangle->GetPointIds()->SetId(1, mesh2D.Cells[i].v1);
      triangle->GetPointIds()->SetId(2, mesh2D.Cells[i].v2);
      cellArray->InsertNextCell(triangle);
    }

    // Set grid data
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_TRIANGLE, cellArray);

    // Write file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
  }

  /// Write 3D surface to VTK (.vtu) file.
  ///
  /// @param surface3D The 3D surface
  /// @parame fileName Filename (path)
  static void Write(const Surface3D &surface3D, std::string fileName)
  {
    info("VTK: Writing 3D surface mesh to file " + fileName);

    // Set vertex data
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (uint i = 0; i < surface3D.Vertices.size(); i++)
    {
      points->InsertNextPoint(surface3D.Vertices[i].x, surface3D.Vertices[i].y,
                              surface3D.Vertices[i].z);
    }

    // Set cell data
    vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
    vtkSmartPointer<vtkCellArray> cellArray =
        vtkSmartPointer<vtkCellArray>::New();
    for (uint i = 0; i < surface3D.Faces.size(); i++)

    {
      triangle->GetPointIds()->SetId(0, surface3D.Faces[i].v0);
      triangle->GetPointIds()->SetId(1, surface3D.Faces[i].v1);
      triangle->GetPointIds()->SetId(2, surface3D.Faces[i].v2);
      cellArray->InsertNextCell(triangle);
    }

    // Set grid data
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_TRIANGLE, cellArray);

    // Write file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
  }

  /// Write 2D gridfield to VTK (.vts) file.
  ///
  /// @param gridfield2D The 2D gridfield
  /// @parame fileName Filename (path)
  static void Write(const GridField2D &gridField2D, std::string fileName)
  {
    info("VTK: Writing 2D grid field to file " + fileName);

    // Set grid data
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkStructuredGrid> structuredGrid =
        vtkSmartPointer<vtkStructuredGrid>::New();
    structuredGrid->SetDimensions(gridField2D.Grid.XSize,
                                  gridField2D.Grid.YSize, 1);
    const Point2D &P = gridField2D.Grid.BoundingBox.P;
    for (uint j = 0; j < gridField2D.Grid.YSize; ++j)
    {
      for (uint i = 0; i < gridField2D.Grid.XSize; ++i)
      {
        points->InsertNextPoint(P.x + i * gridField2D.Grid.XStep,
                                P.y + j * gridField2D.Grid.YStep, 0);
      }
    }
    structuredGrid->SetPoints(points);

    // Set values
    vtkFloatArray *Values = vtkFloatArray::New();
    Values->SetName("Values");
    for (uint i = 0; i < gridField2D.Values.size(); i++)
      Values->InsertTuple1(i, gridField2D.Values[i]);
    structuredGrid->GetPointData()->SetScalars(Values);

    // Write file
    vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
  }

  /// Write 3D gridfield to VTK (.vts) file.
  ///
  /// @param gridfield3D The 3D gridfield
  /// @parame fileName Filename (path)
  static void Write(const GridField3D &gridField3D, std::string fileName)
  {
    info("VTK: Writing 3D grid field to file " + fileName);

    // Set grid data
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkStructuredGrid> structuredGrid =
        vtkSmartPointer<vtkStructuredGrid>::New();
    structuredGrid->SetDimensions(
        gridField3D.Grid.XSize, gridField3D.Grid.YSize, gridField3D.Grid.ZSize);
    const Point3D &P = gridField3D.Grid.BoundingBox.P;
    for (uint k = 0; k < gridField3D.Grid.ZSize; ++k)
    {
      for (uint j = 0; j < gridField3D.Grid.YSize; ++j)
      {
        for (uint i = 0; i < gridField3D.Grid.XSize; ++i)
        {
          points->InsertNextPoint(P.x + i * gridField3D.Grid.XStep,
                                  P.y + j * gridField3D.Grid.YStep,
                                  P.z + k * gridField3D.Grid.ZStep);
        }
      }
    }
    structuredGrid->SetPoints(points);

    // Set values
    vtkFloatArray *Values = vtkFloatArray::New();
    Values->SetName("Values");
    for (uint i = 0; i < gridField3D.Values.size(); i++)
      Values->InsertTuple1(i, gridField3D.Values[i]);
    structuredGrid->GetPointData()->SetScalars(Values);

    // Write file
    vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
  }
};

} // namespace DTCC

#endif
