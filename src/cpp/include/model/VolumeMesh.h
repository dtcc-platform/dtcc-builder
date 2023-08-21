// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_VOLUME_MESH_H
#define DTCC_VOLUME_MESH_H

#include <vector>

#include "Logging.h"
#include "model/Point.h"
#include "model/Simplices.h"
#include "model/Vector.h"

namespace DTCC_BUILDER
{

/// VolumeMesh represents a tetrahedral mesh in 3D

class VolumeMesh : public Printable
{
public:
  /// Array of vertices
  std::vector<Point3D> Vertices{};

  /// Array of cells (tetrahedra)
  std::vector<Simplex3D> Cells{};

  /// Array of cell markers
  std::vector<int> Markers{};

  size_t num_layers{};

  VolumeMesh() = default;
  virtual ~VolumeMesh() {} // make the destructor virtual

  /// Compute of cell
  Point3D MidPoint(size_t cellIndex) const
  {
    Vector3D c;
    c += Vector3D(Vertices[Cells[cellIndex].v0]);
    c += Vector3D(Vertices[Cells[cellIndex].v1]);
    c += Vector3D(Vertices[Cells[cellIndex].v2]);
    c += Vector3D(Vertices[Cells[cellIndex].v3]);
    c /= 4.0;
    return c;
  }

  // Pretty-print
  std::string __str__() const override
  {
    return "VolumeMesh mesh with " + str(Vertices.size()) + " vertices and " +
           str(Cells.size()) + " cells";
  }
};

} // namespace DTCC_BUILDER

#endif
