// VTK I/O
// Anders Logg 2020
// Licensed under the MIT License

#ifndef DTCC_VTK_H
#define DTCC_VTK_H

#include <fstream>
#include <iostream>
#include <vector>

#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "Mesh.h"
#include "Surface.h"
#include "Vector.h"
#include "Simplex.h"
#include "Logging.h"

namespace DTCC
{

class VTK
{
public:
  static void Write(const Mesh3D &mesh3D, std::string fileName)
  {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    for (uint i = 0; i < mesh3D.Vertices.size(); i++)
    {
      points->InsertNextPoint(mesh3D.Vertices[i].x, mesh3D.Vertices[i].y,
                              mesh3D.Vertices[i].z);
    }

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

  // Write 2D mesh to VTK file
  static void Write(const Mesh2D &mesh, const std::string& fileName)
  {
    Info("VTK: Writing 2D mesh to file " + fileName);

    // Open file
    std::ofstream f(fileName.c_str());

    // Check file
    if (!f)
      throw std::runtime_error("Unable to write to file: " + fileName);

    // Set precision
    f << std::setprecision(Parameters::Precision);

    // Write header
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
    f << "<UnstructuredGrid>\n";
    f << "<Piece NumberOfPoints=\"" << mesh.Vertices.size() << "\""
      << " NumberOfCells=\"" << mesh.Cells.size() << "\">" << "\n";

    // Write points
    f << "<Points>\n";
    f << "<DataArray  type=\"Float64\" NumberOfComponents=\"2\" format=\"ascii\">";
    for (const auto& p: mesh.Vertices)
      f << p.x << " " << p.y << " ";
    f << "\n</DataArray>\n";
    f << "</Points>\n";

    // Write cells
    f << "<Cells>\n";
    f << "<DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& c: mesh.Cells)
      f << c.v0 << " " << c.v1 << " " << c.v2 << " ";
    f << "\n</DataArray>\n";
    f << "<DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t i = 0; i < mesh.Cells.size(); i++)
      f << 3*(i + 1) << " ";
    f << "\n</DataArray>\n";
    f << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < mesh.Cells.size(); i++)
      f << "5 ";
    f << "\n</DataArray>\n";
    f << "</Cells>\n";

    // Write footer
    f << "</Piece>\n";
    f << "</UnstructuredGrid>\n";
    f << "</VTKFile>\n";

    // Close file
    f.close();
  }

  // Write 3D mesh to VTK file
  static void Write(const Mesh3D &mesh, const std::string& fileName)
  {
    Info("VTK: Writing 3D mesh to file " + fileName);

    // Open file
    std::ofstream f(fileName.c_str());

    // Check file
    if (!f)
      throw std::runtime_error("Unable to write to file: " + fileName);

    // Set precision
    f << std::setprecision(Parameters::Precision);

    // Write header
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
    f << "<UnstructuredGrid>\n";
    f << "<Piece NumberOfPoints=\"" << mesh.Vertices.size() << "\""
      << " NumberOfCells=\"" << mesh.Cells.size() << "\">\n";

    // Write points
    f << "<Points>\n";
    f << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& p: mesh.Vertices)
      f << p.x << " " << p.y << " " << p.z << " ";
    f << "\n</DataArray>\n";
    f << "</Points>\n";

    // Write cells
    f << "<Cells>\n";
    f << "<DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& c: mesh.Cells)
      f << c.v0 << " " << c.v1 << " " << c.v2 << " " << c.v3 << " ";
    f << "\n</DataArray>\n";
    f << "<DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t i = 0; i < mesh.Cells.size(); i++)
      f << 4*(i + 1) << " ";
    f << "\n</DataArray>\n";
    f << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < mesh.Cells.size(); i++)
      f << "10 ";
    f << "\n</DataArray>\n";
    f << "</Cells>\n";

    // Write footer
    f << "</Piece>\n";
    f << "</UnstructuredGrid>\n";
    f << "</VTKFile>\n";

    // Close file
    f.close();
  }

  // Write 3D surface to VTK file
  static void Write(const Surface3D &surface, const std::string& fileName)
  {
    Info("VTK: Writing 3D surface to file " + fileName);

    // Open file
    std::ofstream f(fileName.c_str());

    // Check file
    if (!f)
      throw std::runtime_error("Unable to write to file: " + fileName);

    // Set precision
    f << std::setprecision(Parameters::Precision);

    // Write header
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
    f << "<UnstructuredGrid>\n";
    f << "<Piece NumberOfPoints=\"" << surface.Vertices.size() << "\""
      << " NumberOfCells=\"" << surface.Cells.size() << "\">\n";

    // Write points
    f << "<Points>\n";
    f << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& p: surface.Vertices)
      f << p.x << " " << p.y << " " << p.z << " ";
    f << "\n</DataArray>\n";
    f << "</Points>\n";

    // Write cells
    f << "<Cells>\n";
    f << "<DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& c: surface.Cells)
      f << c.v0 << " " << c.v1 << " " << c.v2 << " ";
    f << "\n</DataArray>\n";
    f << "<DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t i = 0; i < surface.Cells.size(); i++)
      f << 3*(i + 1) << " ";
    f << "\n</DataArray>\n";
    f << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < surface.Cells.size(); i++)
      f << "5 ";
    f << "\n</DataArray>\n";
    f << "</Cells>\n";

    // Write footer
    f << "</Piece>\n";
    f << "</UnstructuredGrid>\n";
    f << "</VTKFile>\n";

    // Close file
    f.close();
  }

};

} // namespace DTCC

#endif
