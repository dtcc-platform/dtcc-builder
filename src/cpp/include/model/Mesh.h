// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_H
#define DTCC_MESH_H

#include <vector>

#include "Logging.h"
#include "model/Point.h"
#include "model/Simplices.h"
#include "model/Vector.h"

namespace DTCC_BUILDER
{

/// Mesh represents a triangular mesh in 3D

class Mesh : public Printable
{
public:
  /// Array of vertices
  std::vector<Point3D> Vertices{};

  /// Array of faces (triangles)
  std::vector<Simplex2D> Faces{};

  /// Array of normals
  std::vector<Point3D> Normals{};

  /// Array of cell markers
  std::vector<int> Markers{};

  Mesh() = default;
  virtual ~Mesh() {} // make the destructor virtual

  /// Compute midpoint of cell
  Point3D MidPoint(size_t cellIndex) const
  {
    Vector3D c{};
    c += Vector3D(Vertices[Faces[cellIndex].v0]);
    c += Vector3D(Vertices[Faces[cellIndex].v1]);
    c += Vector3D(Vertices[Faces[cellIndex].v2]);
    c /= 3.0;
    return c;
  }

  /// Pretty-print
  std::string __str__() const override
  {
    return "Mesh with " + str(Vertices.size()) + " vertices and " +
           str(Faces.size()) + " faces";
  }
};

} // namespace DTCC_BUILDER

#endif
