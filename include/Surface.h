// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_SURFACE_H
#define DTCC_SURFACE_H

#include <vector>

#include "Logging.h"
#include "Point.h"
#include "Simplex.h"

namespace DTCC
{

  class Surface2D : public Printable
  {
  public:

    // Array of vertices
    std::vector<Point2D> Vertices;

    // Array of edges
    std::vector<Simplex1D> Edges;

    // Array of normals
    std::vector<Vector2D> Normals;

    /// Pretty-print
    std::string __str__() const override
    {
      return "2D surface (mesh boundary) with " + str(Vertices.size()) +
             " vertices and " + str(Edges.size()) + " edges";
    }

  };

  class Surface3D : public Printable
  {
  public:

    // Array of vertices
    std::vector<Point3D> Vertices;

    // Array of faces
    std::vector<Simplex2D> Faces;

    // Array of normal
    std::vector<Vector3D> Normals;

    /// Pretty-print
    std::string __str__() const override
    {
      return "3D surface (mesh boundary) with " + str(Vertices.size()) +
             " vertices and " + str(Faces.size()) + " faces";
    }

  };

} // namespace DTCC

#endif
