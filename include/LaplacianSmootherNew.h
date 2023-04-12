// Copyright (C) 2023 Authors
// Licensed under the MIT License

#ifndef DTCC_LAPLACIAN_SMOOTHER_NEW_H
#define DTCC_LAPLACIAN_SMOOTHER_NEW_H

#include "Mesh.h"
#include "Timer.h"

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

    // FIXME: Write code here
  }
};

} // namespace DTCC

#endif // DTCC_LAPLACIAN_SMOOTHER_NEW_H
