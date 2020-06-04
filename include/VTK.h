// VTK I/O
// Anders Logg 2020

#ifndef DTCC_VTK_H
#define DTCC_VTK_H

#include <fstream>
#include <iostream>
#include <vector>

#include "Mesh.h"
#include "Surface.h"
#include "Point.h"
#include "Simplex.h"

namespace DTCC
{

class VTK
{
public:

  // Write 2D mesh to VTK file
  static void Write(const Mesh2D &mesh, std::string fileName)
  {
    std::cout << "VTK: Writing 2D mesh to file " << fileName << std::endl;

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
  static void Write(const Mesh3D &mesh, std::string fileName)
  {
    std::cout << "VTK: Writing 3D mesh to file " << fileName << std::endl;

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
  static void Write(const Surface3D &surface, std::string fileName)
  {
    std::cout << "VTK: Writing 3D surface to file " << fileName << std::endl;

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
