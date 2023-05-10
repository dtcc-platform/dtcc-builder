// Copyright (C) 2023 Authors
// Licensed under the MIT License

// cd dtcc-builder/build
// make dtcc-generate-mesh
// dtcc-generate-mesh/dtcc-generate-mesh
// ../data/HelsingborgResidential2022/parameters-new.json
//
// New Laplacian Smoother class
// TO DO:
//
// 1) Fix wavey sidewalls bug.
//    Probably caused by bad BCs

#ifndef DTCC_LAPLACIAN_SMOOTHER_NEW_H
#define DTCC_LAPLACIAN_SMOOTHER_NEW_H

#include "Mesh.h"
#include "Timer.h"

#include "../sandbox/smoothing-2023/include/boundaryConditions.hpp"
#include "../sandbox/smoothing-2023/include/stiffnessMatrix.hpp"

#define MAX_ITER 5000

namespace DTCC
{

bool checkSidewallVertices(const Point3D &vertex, const int Markers);

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
    info(mesh3D.__str__());
    Timer timer("SmoothMesh3DNew");

    // Local Stifness Matrices
    stiffnessMatrix AK(mesh3D);

    // Solution Vector
    std::vector<double> u(mesh3D.Vertices.size(), 0);
    // Load Vector
    std::vector<double> b(mesh3D.Vertices.size(), 0);

    BoundaryConditions bc(mesh3D, cityModel, dem, topHeight, fixBuildings);
    bc.apply(b);
    bc.apply(AK);

    // Initial Approximation of the solution
    if (!fixBuildings)
    {
      u = initialGuess(mesh3D, dem, topHeight, bc);
    }
    else
    {
      u = b;
    }

    // UnassembledJacobi(mesh3D, AK, b, u);
    UnassembledGaussSeidel(mesh3D, AK, b, u);

    double minElevation = dem.Min();
    size_t c1 = 0;
    size_t c2 = 0;

    // Update mesh coordinates
    for (std::size_t i = 0; i < mesh3D.Vertices.size(); i++)
    {
      mesh3D.Vertices[i].z += u[i];
      if (mesh3D.Vertices[i].z > topHeight)
      {
        std::cout << i << "z: " << mesh3D.Vertices[i].z << std::endl;
        c1++;
      }
      if (mesh3D.Vertices[i].z < minElevation)
      {
        std::cout << i << "z: " << mesh3D.Vertices[i].z << " " << u[i]
                  << std::endl;
        c2++;
      }
    }
    std::cout << "Over TOP: c1 = " << c1 << "\nBelow Min: c2 = " << c2
              << " MIN: " << minElevation << " MAX " << dem.Max() << std::endl;
  }

  static void UnassembledJacobi(const Mesh3D &mesh3D,
                                stiffnessMatrix &A,
                                std::vector<double> &b,
                                std::vector<double> &u,
                                const size_t max_iterations = MAX_ITER,
                                const double tolerance = 1e-16)
  {
    info("Element-by-Element Jacobi Solver");
    Timer timer("EbE Jacobi");

    const size_t nV = mesh3D.Vertices.size();
    const size_t nC = mesh3D.Cells.size();

    // Non-diagonal Elements sum
    std::vector<double> c(nV);

    std::array<uint, 4> I = {0};

    size_t iterations;
    double residual;
    for (iterations = 0; iterations < max_iterations; iterations++)
    {
      c = b;
      residual = 0;
      for (size_t cn = 0; cn < nC; cn++)
      {
        I[0] = mesh3D.Cells[cn].v0;
        I[1] = mesh3D.Cells[cn].v1;
        I[2] = mesh3D.Cells[cn].v2;
        I[3] = mesh3D.Cells[cn].v3;
        for (u_int8_t i = 0; i < 4; i++)
        {
          c[I[i]] -= A(cn, i, (i + 1) % 4) * u[I[(i + 1) % 4]] +
                     A(cn, i, (i + 2) % 4) * u[I[(i + 2) % 4]] +
                     A(cn, i, (i + 3) % 4) * u[I[(i + 3) % 4]];
        }
      }

      // Update solution vector
      for (size_t i = 0; i < nV; i++)
      {
        double res = u[i];
        u[i] = c[i] / A.diagonal[i];
        res = std::abs(res - u[i]);
        residual = std::max(residual, res);
      }

      // Check Convergance
      if (residual < tolerance)
        break;
    }
    timer.Stop();
    timer.Print();

    std::cout << "Jacobi finished after " << iterations << " / "
              << max_iterations << " iterations" << std::endl;
    std::cout << "With residual: " << residual << std::endl;
    // std::cout << "Execution Time:" << ms_double.count() << "ms " <<
    // std::endl;
  }

  // Compute the number of Cells that each Vertex belongs
  static void _getVertexDegree(std::vector<uint> &VertexDegree,
                               const Mesh3D &mesh)
  {
    for (size_t cn = 0; cn < mesh.Cells.size(); cn++)
    {
      VertexDegree[mesh.Cells[cn].v0]++;
      VertexDegree[mesh.Cells[cn].v1]++;
      VertexDegree[mesh.Cells[cn].v2]++;
      VertexDegree[mesh.Cells[cn].v3]++;
    }
  }

  static void UnassembledGaussSeidel(const Mesh3D &mesh3D,
                                     stiffnessMatrix &A,
                                     std::vector<double> &b,
                                     std::vector<double> &u,
                                     const size_t max_iterations = MAX_ITER,
                                     const double tolerance = 1e-16)
  {
    info("Element-by-Element Gauss-Seidel Solver");
    Timer timer("EbE Gauss Seidel");

    const size_t nV = mesh3D.Vertices.size();
    const size_t nC = mesh3D.Cells.size();

    // Non-diagonal Elements sum
    std::vector<double> c(nV);

    std::array<uint, 4> I = {0};

    // Compute the number of Cells that each Vertex belongs
    std::vector<uint> vertexDegree(nV);
    std::vector<uint> vDeg(nV);
    _getVertexDegree(vertexDegree, mesh3D);

    size_t iterations;
    double residual;
    for (iterations = 0; iterations < max_iterations; iterations++)
    {
      c = b;
      vDeg = vertexDegree;
      residual = 0;
      for (size_t cn = 0; cn < nC; cn++)
      {
        I[0] = mesh3D.Cells[cn].v0;
        I[1] = mesh3D.Cells[cn].v1;
        I[2] = mesh3D.Cells[cn].v2;
        I[3] = mesh3D.Cells[cn].v3;
        for (u_int8_t i = 0; i < 4; i++)
        {
          c[I[i]] -=
              A._data[cn * 16 + i * 4 + (i + 1) % 4] * u[I[(i + 1) % 4]] +
              A._data[cn * 16 + i * 4 + (i + 2) % 4] * u[I[(i + 2) % 4]] +
              A._data[cn * 16 + i * 4 + (i + 3) % 4] * u[I[(i + 3) % 4]];

          // c[I[i]] -= A(cn,i,(i + 1) % 4) * u[I[(i + 1) % 4]] +
          //            A(cn,i,(i + 2) % 4) * u[I[(i + 2) % 4]] +
          //            A(cn,i,(i + 3) % 4) * u[I[(i + 3) % 4]] ;

          vDeg[I[i]]--;
          if (vDeg[I[i]] == 0)
          {
            double res = u[I[i]];
            u[I[i]] = c[I[i]] / A.diagonal[I[i]];
            res = std::abs(res - u[I[i]]);
            residual = std::max(residual, res);
          }
        }
      }
      // std::cout << iterations << "  err: " << error << std::endl;

      // Check Convergance
      if (residual < tolerance)
        break;
    }
    timer.Stop();
    timer.Print();

    std::cout << "Gauss-Seidel finished after " << iterations << " / "
              << max_iterations << " iterations" << std::endl;
    std::cout << "With residual: " << residual << std::endl;
    // std::cout << "Execution Time:" << ms_double.count() << "ms " <<
    // std::endl;
  }

  // Placeholder Initial Guess
  static std::vector<double> initialGuess(const Mesh3D &mesh3D,
                                          const GridField2D &dem,
                                          double topHeight,
                                          BoundaryConditions &bc)
  {
    info("Initial Approximation for Solution Vector");
    const std::size_t nV = mesh3D.Vertices.size();

    std::vector<double> u = bc.values;

    double meanElevation = dem.Mean();
    double maxElevation = dem.Max();
    for (size_t i = 0; i < nV; i++)
    {
      if (bc.vMarkers[i] == -4)
      {
        // if(checkSidewallVertices(mesh3D.Vertices[i],bc.vMarkers[i]))
        // {
        //   u[i] =meanElevation * (1 - mesh3D.Vertices[i].z / topHeight);
        // }else
        {
          const Vector2D p(mesh3D.Vertices[i].x, mesh3D.Vertices[i].y);
          u[i] = dem(p) * (1 - mesh3D.Vertices[i].z / topHeight);
        }
      }
    }
    return u;
  }
};

// Test:: Probably Gonna Remove this!
bool checkSidewallVertices(const Point3D &vertex, int Marker)
{
  const double x_min = 0.0;
  const double x_max = 999.99;

  const double y_min = 0.00255543;
  const double y_max = 995.833;

  const double epsilon = 0.01;

  bool onDomainboundary = false;

  if (abs(vertex.x - x_min) < epsilon)
  {
    onDomainboundary = true;
    // std::cout << "Boundary Vertex: " << " x: " << vertex.x
    //           << " y: " << vertex.y << " z : " << vertex.z << " Marker: " <<
    //           Marker << std::endl;
  }
  else if (abs(vertex.x - x_max) < epsilon)
  {
    onDomainboundary = true;
    // std::cout << " Boundary Vertex: " << " x: " << vertex.x
    //           << " y: " << vertex.y << " z : " << vertex.z << " Marker: " <<
    //           Marker
    //           << std::endl;
  }
  else if (abs(vertex.y - y_min) < epsilon)
  {
    onDomainboundary = true;
    // std::cout << "Boundary Vertex: " << " x: " << vertex.x
    //           << " y: " << vertex.y << " z : " << vertex.z << " Marker: " <<
    //           Marker
    //           << std::endl;
  }
  else if (abs(vertex.y - y_max) < epsilon)
  {
    onDomainboundary = true;
    // std::cout << "Boundary Vertex: " << " x: " << vertex.x
    //           << " y: " << vertex.y << " z : " << vertex.z << " Marker: " <<
    //           Marker
    //           << std::endl;
  }
  return onDomainboundary;
}

} // namespace DTCC

// namespace DTCC

#endif // DTCC_LAPLACIAN_SMOOTHER_NEW_H
