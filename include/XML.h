// FEniCS XML I/O
// Anders Logg 2019
// Licensed under the MIT License

#ifndef DTCC_XML_H
#define DTCC_XML_H

#include <fstream>
#include <iostream>
#include <vector>

#include "Mesh.h"
#include "Vector.h"
#include "Simplex.h"

namespace DTCCBUILDER
{

class XML
{
public:
  // Write 2D mesh to FEniCS XML file (.xml will be appended)
  static void Write(const Mesh2D &Mesh, std::string Prefix)
  {
    // Open file
    std::string FileName = Prefix + ".xml";
    std::ofstream f(FileName.c_str());

    // Check file
    if (!f)
      throw std::runtime_error("Unable to write to file: " + FileName);

    // Write header
    f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    f << "<dolfin xmlns:dolfin=\"http://www.fenics.org/dolfin/\">\n";
    f << "  <mesh celltype=\"triangle\" dim=\"2\">\n";

    // Write points
    f << "    <vertices size=\"" << Mesh.Points.size() << "\">\n";
    for (size_t i = 0; i < Mesh.Points.size(); i++)
      Write(Mesh.Points[i], i, f);
    f << "    </vertices>\n";

    // Write cells
    f << "    <cells size=\"" << Mesh.Cells.size() << "\">\n";
    for (size_t i = 0; i < Mesh.Cells.size(); i++)
      Write(Mesh.Cells[i], i, f);
    f << "    </cells>\n";

    // Write footer
    f << "  </mesh>\n";
    f << "</dolfin>\n";

    // Close file
    f.close();
  }

  // Write 3D mesh to FEniCS XML file (.xml will be appended)
  static void Write(const Mesh3D &Mesh, std::string Prefix)
  {
    // Open file
    std::string FileName = Prefix + ".xml";
    std::ofstream f(FileName.c_str());

    // Check file
    if (!f)
      throw std::runtime_error("Unable to write to file: " + FileName);

    // Set precision
    f << std::setprecision(Parameters::Precision);

    // Write header
    f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    f << "<dolfin xmlns:dolfin=\"http://www.fenics.org/dolfin/\">\n";
    f << "  <mesh celltype=\"tetrahedron\" dim=\"3\">\n";

    // Write points
    f << "    <vertices size=\"" << Mesh.Points.size() << "\">\n";
    for (size_t i = 0; i < Mesh.Points.size(); i++)
      Write(Mesh.Points[i], i, f);
    f << "    </vertices>\n";

    // Write cells
    f << "    <cells size=\"" << Mesh.Cells.size() << "\">\n";
    for (size_t i = 0; i < Mesh.Cells.size(); i++)
      Write(Mesh.Cells[i], i, f);
    f << "    </cells>\n";

    // Write footer
    f << "  </mesh>\n";
    f << "</dolfin>\n";

    // Close file
    f.close();
  }

private:
  // Write 2D point to file
  static void Write(const Vector2D &p, size_t i, std::ofstream &f)
  {
    f << "      <vertex"
      << " index=\"" << i << "\""
      << " x=\"" << p.x << "\""
      << " y=\"" << p.y << "\""
      << "/>\n";
  }

  // Write 3D point to file
  static void Write(const Vector3D &p, size_t i, std::ofstream &f)
  {
    f << "      <vertex"
      << " index=\"" << i << "\""
      << " x=\"" << p.x << "\""
      << " y=\"" << p.y << "\""
      << " z=\"" << p.z << "\""
      << "/>\n";
  }

  // Write 2D simplex to file
  static void Write(const Simplex2D &t, size_t i, std::ofstream &f)
  {
    f << "      <triangle"
      << " index=\"" << i << "\""
      << " v0=\"" << t.v0 << "\""
      << " v1=\"" << t.v1 << "\""
      << " v2=\"" << t.v2 << "\""
      << "/>\n";
  }

  // Write 3D simplex to file
  static void Write(const Simplex3D &t, size_t i, std::ofstream &f)
  {
    f << "      <tetrahedron"
      << " index=\"" << i << "\""
      << " v0=\"" << t.v0 << "\""
      << " v1=\"" << t.v1 << "\""
      << " v2=\"" << t.v2 << "\""
      << " v3=\"" << t.v3 << "\""
      << "/>\n";
  }
};

} // namespace DTCCBUILDER

#endif