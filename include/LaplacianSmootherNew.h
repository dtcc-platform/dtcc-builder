// Copyright (C) 2023 Authors
// Licensed under the MIT License

#ifndef DTCC_LAPLACIAN_SMOOTHER_NEW_H
#define DTCC_LAPLACIAN_SMOOTHER_NEW_H

#include "Mesh.h"
#include "Timer.h"

#include "../sandbox/smoothing-2023/include/boundaryConditions.hpp"

namespace DTCC
{

class LaplacianSmootherNew
{
public:
  // Smooth mesh using Laplacian smoothing
  static void SmoothMesh3D(Mesh3D &mesh3D,
                           const CityModel &cityModel,
                           const GridField2D &dem,
                           double topHeight,
                           bool fixBuildings)
  {
    info("LaplacianSmoother: Smoothing mesh (Laplacian smoothing NEW)...");
    Timer timer("SmoothMesh3DNew");

    std::cout << mesh3D.Markers.size() << std::endl;

    int *vMarkers = new int[mesh3D.Vertices.size()];
    getVerticeMarkers(mesh3D, vMarkers);
  }
};

} // namespace DTCC

#endif // DTCC_LAPLACIAN_SMOOTHER_NEW_H
