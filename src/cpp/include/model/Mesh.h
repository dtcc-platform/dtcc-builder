// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_H
#define DTCC_MESH_H

#include <vector>

#include "Logging.h"
#include "model/Simplices.h"
#include "model/Vector.h"

namespace DTCC_BUILDER
{

/// Mesh represents a triangular mesh in 3D

class Mesh : public Printable
{
public:
  /// Array of vertices
  std::vector<Vector3D> vertices{};

  /// Array of faces (triangles)
  std::vector<Simplex2D> faces{};

  /// Array of normals
  std::vector<Vector3D> normals{};

  /// Array of cell markers
  std::vector<int> markers{};

  Mesh() = default;
  virtual ~Mesh() {} // make the destructor virtual

  /// Compute midpoint of cell
  Vector3D mid_point(size_t cell_index) const
  {
    Vector3D c{};
    c += Vector3D(vertices[faces[cell_index].v0]);
    c += Vector3D(vertices[faces[cell_index].v1]);
    c += Vector3D(vertices[faces[cell_index].v2]);
    c /= 3.0;
    return c;
  }

  /// Pretty-print
  std::string __str__() const override
  {
    return "Mesh with " + str(vertices.size()) + " vertices and " +
           str(faces.size()) + " faces";
  }
};

} // namespace DTCC_BUILDER

#endif
