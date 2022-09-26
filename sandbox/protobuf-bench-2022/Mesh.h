// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_H
#define DTCC_MESH_H

#include <vector>

#include "Logging.h"
#include "Point.h"
#include "Simplex.h"
#include "Vector.h"

namespace DTCC
{


  class Mesh3D
  {
  public:

    /// Array of vertices
    std::vector<Point3D> Vertices{};

    /// Array of cells (tetrahedra)
    std::vector<Simplex3D> Cells{};

    /// Array of cell markers
    std::vector<int> Markers{};

  };

} // namespace DTCC

#endif
