// Copyright (C) 2023 George Spaias
// Licensed under the MIT License

#ifndef DTCC_SMOOTHER_H
#define DTCC_SMOOTHER_H

#include "BoundaryConditions.h"
#include "Mesh.h"
#include "StiffnessMatrix.h"
#include "Timer.h"

namespace DTCC_BUILDER
{

bool checkSidewallVertices(const Point3D &vertex, const int Markers);

class Smoother
{
public:
  // Smooth mesh using Laplacian smoothing
  static void smooth_volume_mesh(Mesh3D &volume_mesh,
                                 const CityModel &city_model,
                                 const GridField2D &dem,
                                 double top_height,
                                 bool fix_buildings,
                                 size_t max_iterations,
                                 double relative_tolerance)
  {
    info("LaplacianSmoother: Smoothing mesh (Laplacian smoothing NEW)...");
    info(volume_mesh.__str__());

    // Local Stifness Matrices
    stiffnessMatrix AK(volume_mesh);

    // Solution vector and load vector
    std::vector<double> u(volume_mesh.Vertices.size(), 0);
    std::vector<double> b(volume_mesh.Vertices.size(), 0);

    BoundaryConditions bc(volume_mesh, city_model, dem, top_height,
                          fix_buildings);
    bc.apply(b);
    bc.apply(AK);

    // Initial Approximation of the solution
    if (!fix_buildings)
      u = initialGuess(volume_mesh, dem, top_height, bc);
    else
      u = b;

    // UnassembledJacobi(mesh3D, AK, b, u);
    UnassembledGaussSeidel(volume_mesh, AK, b, u, max_iterations,
                           relative_tolerance);

    double minElevation = dem.Min();
    size_t c1 = 0;
    size_t c2 = 0;

    // Update mesh coordinates
    for (std::size_t i = 0; i < volume_mesh.Vertices.size(); i++)
      volume_mesh.Vertices[i].z += u[i];
  }

  static void UnassembledJacobi(const Mesh3D &mesh3D,
                                stiffnessMatrix &A,
                                std::vector<double> &b,
                                std::vector<double> &u,
                                const size_t maxIter = 1000,
                                const double relTol = 1e-16)
  {
    info("Element-by-Element Jacobi Solver");
    // Timer timer("EbE Jacobi");

    const size_t nV = mesh3D.Vertices.size();
    const size_t nC = mesh3D.Cells.size();

    // Non-diagonal Elements sum
    std::vector<double> c(nV);

    std::array<uint, 4> I = {0};

    size_t iterations;
    double residual;
    for (iterations = 0; iterations < maxIter; iterations++)
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
      if (residual < relTol)
        break;
    }

    std::cout << "Jacobi finished after " << iterations << " / " << maxIter
              << " iterations" << std::endl;
    std::cout << "With residual: " << residual << std::endl;
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
                                     const size_t maxIter,
                                     const double relTol = 1e-16)
  {
    info("Element-by-Element Gauss-Seidel Solver");

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
    for (iterations = 0; iterations < maxIter; iterations++)
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

      // Check Convergance
      if (residual < relTol)
        break;
    }

    std::cout << "Gauss-Seidel finished after " << iterations << " / "
              << maxIter << " iterations" << std::endl;
    std::cout << "With residual: " << residual << std::endl;
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
        const Vector2D p(mesh3D.Vertices[i].x, mesh3D.Vertices[i].y);
        u[i] = dem(p) * (1 - mesh3D.Vertices[i].z / topHeight);
      }
    }
    return u;
  }
};

} // namespace DTCC_BUILDER

#endif // DTCC_LAPLACIAN_SMOOTHER_NEW_H
