// FEniCS utiliy functions.
// Copyright (C) 2019 Anders Logg.

#ifndef VC_FENICS_H
#define VC_FENICS_H

#include <dolfin.h>
#include "Mesh.h"

namespace VirtualCity
{

class FEniCS
{
public:
  // Create FEniCS mesh from VirtualCity mesh (2D)
  static void ConvertMesh(const Mesh2D &mesh2D, dolfin::Mesh &mesh)
  {
    // Inialize mesh editor
    dolfin::MeshEditor meshEditor;
    meshEditor.open(mesh, "triangle", 2, 2);

    // Add vertices
    meshEditor.init_vertices(mesh2D.Points.size());
    for (size_t i = 0; i < mesh2D.Points.size(); i++)
    {
      const Point2D &p = mesh2D.Points[i];
      meshEditor.add_vertex(i, p.x, p.y);
    }

    // Add cells
    meshEditor.init_cells(mesh2D.Cells.size());
    for (size_t i = 0; i < mesh2D.Cells.size(); i++)
    {
      const Simplex2D &c = mesh2D.Cells[i];
      meshEditor.add_cell(i, c.v0, c.v1, c.v2);
    }

    // Finalize mesh editor
    meshEditor.close();
  }

  // Create FEniCS mesh from VirtualCity mesh (3D)
  static void ConvertMesh(const Mesh3D &mesh3D, dolfin::Mesh &mesh)
  {
    // Inialize mesh editor
    dolfin::MeshEditor meshEditor;
    meshEditor.open(mesh, "tetrahedron", 3, 3);

    // Add vertices
    meshEditor.init_vertices(mesh3D.Points.size());
    for (size_t i = 0; i < mesh3D.Points.size(); i++)
    {
      const Point3D &p = mesh3D.Points[i];
      meshEditor.add_vertex(i, p.x, p.y, p.z);
    }

    // Add cells
    meshEditor.init_cells(mesh3D.Cells.size());
    for (size_t i = 0; i < mesh3D.Cells.size(); i++)
    {
      const Simplex3D &c = mesh3D.Cells[i];
      meshEditor.add_cell(i, c.v0, c.v1, c.v2, c.v3);
    }

    // Finalize mesh editor
    meshEditor.close();
  }
};

} // namespace VirtualCity

#endif
