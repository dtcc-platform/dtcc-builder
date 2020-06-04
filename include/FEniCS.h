// FEniCS utiliy functions.
// Copyright (C) 2019 Anders Logg.

#ifndef DTCC_FENICS_H
#define DTCC_FENICS_H

#include "Mesh.h"
#include "Surface.h"
#include <dolfin.h>

namespace DTCC
{

class FEniCS
{
public:
  // Create FEniCS mesh from DTCC mesh (2D)
  static void ConvertMesh(const Mesh2D &mesh2D, dolfin::Mesh &mesh)
  {
    // Inialize mesh editor
    dolfin::MeshEditor meshEditor;
    meshEditor.open(mesh, "triangle", 2, 2);

    // Add vertices
    meshEditor.init_vertices(mesh2D.Vertices.size());
    for (size_t i = 0; i < mesh2D.Vertices.size(); i++)
    {
      const Point2D &p = mesh2D.Vertices[i];
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

  // Create FEniCS mesh from DTCC mesh (3D)
  static void ConvertMesh(const Mesh3D &mesh3D, dolfin::Mesh &mesh)
  {
    // Inialize mesh editor
    dolfin::MeshEditor meshEditor;
    meshEditor.open(mesh, "tetrahedron", 3, 3);

    // Add vertices
    meshEditor.init_vertices(mesh3D.Vertices.size());
    for (size_t i = 0; i < mesh3D.Vertices.size(); i++)
    {
      const Point3D &p = mesh3D.Vertices[i];
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

  // Create FEniCS mesh from DTCC surface (3D)
  static void ConvertMesh(const Surface3D &surface,
                          dolfin::Mesh &mesh)
  {
    std::vector<Surface3D> surfaces;
    surfaces.push_back(surface);
    ConvertMesh(surfaces, mesh);
  }

  // Create FEniCS mesh from DTCC surfaces (3D)
  static void ConvertMesh(const std::vector<Surface3D> &surfaces,
                          dolfin::Mesh &mesh)
  {
    // Inialize mesh editor
    dolfin::MeshEditor meshEditor;
    meshEditor.open(mesh, "triangle", 2, 3);

    // Count the number of vertices and triangles
    size_t numVertices = 0;
    size_t numTriangles = 0;
    for (auto const &surface : surfaces)
    {
      numVertices += surface.Vertices.size();
      numTriangles += surface.Cells.size();
    }

    // Add vertices
    meshEditor.init_vertices(numVertices);
    {
      size_t k = 0;
      for (auto const &surface : surfaces)
      {
        for (size_t i = 0; i < surface.Vertices.size(); i++)
        {
          const Point3D &p = surface.Vertices[i];
          meshEditor.add_vertex(k++, p.x, p.y, p.z);
        }
      }
    }

    // Add cells
    meshEditor.init_cells(numTriangles);
    {
      size_t k = 0;
      size_t offset = 0;
      for (auto const &surface : surfaces)
      {
        for (size_t i = 0; i < surface.Cells.size(); i++)
        {
          const Simplex2D &c = surface.Cells[i];
          meshEditor.add_cell(k++, c.v0 + offset, c.v1 + offset, c.v2 + offset);
        }
        offset += surface.Vertices.size();
      }
    }

    // Finalize mesh editor
    meshEditor.close();
  }
};

} // namespace DTCC

#endif
