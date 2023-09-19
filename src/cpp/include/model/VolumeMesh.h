// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_VOLUME_MESH_H
#define DTCC_VOLUME_MESH_H

#include <vector>

#include "Logging.h"
#include "model/Simplices.h"
#include "model/Vector.h"

namespace DTCC_BUILDER
{

/// VolumeMesh represents a tetrahedral mesh in 3D

class VolumeMesh : public Printable
{
public:
  /// Array of vertices
  std::vector<Vector3D> vertices{};

  /// Array of cells (tetrahedra)
  std::vector<Simplex3D> cells{};

  /// Array of cell markers
  std::vector<int> markers{};

  size_t num_layers{};

  VolumeMesh() = default;
  virtual ~VolumeMesh() {} // make the destructor virtual

  /// Compute of cell
  Vector3D mid_point(size_t cell_index) const
  {
    Vector3D c;
    c += Vector3D(vertices[cells[cell_index].v0]);
    c += Vector3D(vertices[cells[cell_index].v1]);
    c += Vector3D(vertices[cells[cell_index].v2]);
    c += Vector3D(vertices[cells[cell_index].v3]);
    c /= 4.0;
    return c;
  }

  // Pretty-print
  std::string __str__() const override
  {
    return "VolumeMesh mesh with " + str(vertices.size()) + " vertices and " +
           str(cells.size()) + " cells";
  }
};

} // namespace DTCC_BUILDER

#endif
