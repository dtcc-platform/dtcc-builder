// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_FENICS_H
#define DTCC_FENICS_H

#include <dolfin.h>

#include "Logging.h"
#include "Mesh.h"
#include "Surface.h"

namespace DTCCBUILDER
{

class FEniCS
{
public:
  // Create FEniCS mesh from DTCC mesh (2D)
  static void ConvertMesh(const Mesh2D &mesh2D, dolfin::Mesh &mesh)
  {
    info("FEniCS: Converting to FEniCS 2D mesh...");

    // Inialize mesh editor
    dolfin::MeshEditor meshEditor;
    meshEditor.open(mesh, "triangle", 2, 2);

    // Add vertices
    meshEditor.init_vertices(mesh2D.Vertices.size());
    for (size_t i = 0; i < mesh2D.Vertices.size(); i++)
    {
      const Vector2D &p = Vector2D(mesh2D.Vertices[i]);
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
    info("FEniCS: Converting to FEniCS 3D mesh...");

    // Inialize mesh editor
    dolfin::MeshEditor meshEditor;
    meshEditor.open(mesh, "tetrahedron", 3, 3);

    // Add vertices
    meshEditor.init_vertices(mesh3D.Vertices.size());
    for (size_t i = 0; i < mesh3D.Vertices.size(); i++)
    {
      const Vector3D &p = Vector3D(mesh3D.Vertices[i]);
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
    info("FEniCS: Converting to FEniCS 3D surface...");

    // Inialize mesh editor
    dolfin::MeshEditor meshEditor;
    meshEditor.open(mesh, "triangle", 2, 3);

    // Count the number of vertices and triangles
    size_t numVertices = 0;
    size_t numTriangles = 0;
    for (auto const &surface : surfaces)
    {
      numVertices += surface.Vertices.size();
      numTriangles += surface.Faces.size();
    }

    // Add vertices
    meshEditor.init_vertices(numVertices);
    {
      size_t k = 0;
      for (auto const &surface : surfaces)
      {
        for (size_t i = 0; i < surface.Vertices.size(); i++)
        {
          const Vector3D &p = Vector3D(surface.Vertices[i]);
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
        for (size_t i = 0; i < surface.Faces.size(); i++)
        {
          const Simplex2D &c = surface.Faces[i];
          meshEditor.add_cell(k++, c.v0 + offset, c.v1 + offset, c.v2 + offset);
        }
        offset += surface.Vertices.size();
      }
    }

    // Finalize mesh editor
    meshEditor.close();
  }

  // Write FEniCS mesh to file
  static void Write(const dolfin::Mesh &mesh, const std::string& fileName)
  {
    info("FEniCS: Writing to file " + fileName + "...");
    dolfin::File(fileName) << mesh;
  }

  // Write FEniCS boundary mesh to file
  static void Write(const dolfin::BoundaryMesh &boundary, const std::string& fileName)
  {
    info("FEniCS: Writing to file " + fileName + "...");
    dolfin::File(fileName) << boundary;
  }

  // Write FEniCS function to file
  static void Write(const dolfin::Function &function, const std::string& fileName)
  {
    info("FEniCS: Writing to file " + fileName + "...");
    dolfin::File(fileName) << function;
  }
};

} // namespace DTCCBUILDER

#endif
