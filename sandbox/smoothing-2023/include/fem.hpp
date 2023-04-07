
#ifndef FEM_H
#define FEM_H

#include <algorithm>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "Mesh.h"
#include "datamodel/CityModel.h"

#include "assembled.hpp"
#include "boundaryConditions.hpp"
#include "error.hpp"
#include "stiffnessMatrix.hpp"

#define MAX_ITERATIONS 1000

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

// Test Function| Needs better structure
void smoothLaplaceJacobi(const Mesh3D &mesh,
                         const CityModel &citymodel,
                         const GridField2D &dtm,
                         const size_t max_iterations)
{
  info("Mesh Smoothing with Element-by-Element Jacobi Method");

  // const size_t max_iterations = MAX_ITERATIONS;
  const size_t nC = mesh.Cells.size();
  const size_t nV = mesh.Vertices.size();

  // Compute local Stiffness Matrices
  double *AK = new double[16 * nC];
  compute_transformation_matrix(AK, mesh);

  // Check Boundary Vertices
  bool *isBoundary = new bool[nV]();
  double *b = new double[nV]();
  int *vMarkers = new int[nV];

  // getVerticeMarkers(mesh,vMarkers);

  double topHeight = dtm.Mean() + static_cast<double>(100.0);
  // applyBC_load(vMarkers,mesh, citymodel, dtm,topHeight, b);
  checkBoundaryPoints(isBoundary, b, mesh);

  // Diagonal d
  double *d = new double[nV]();
  // applyBC_stiffnessMat(vMarkers,mesh,AK,d);

  // Initial guess for solution 0 TBD
  double *u = new double[nV]();
  std::memcpy(u, b, nV * sizeof(double));
  // std::memset(u,0,nV * sizeof(double));
  // Non-diagonal Elements sum
  double *c = new double[nV];

  std::array<uint, 4> I = {0};

  for (size_t cn = 0; cn < nC; cn++)
  {
    // Initializing Global Index for each cell
    I[0] = mesh.Cells[cn].v0;
    I[1] = mesh.Cells[cn].v1;
    I[2] = mesh.Cells[cn].v2;
    I[3] = mesh.Cells[cn].v3;

    for (uint8_t i = 0; i < 4; i++)
    {
      if (isBoundary[I[i]] == true)
      {
        d[I[i]] = 1.0;
        memset(&AK[cn * 16 + i * 4], 0, 4 * sizeof(double));
        // for (uint8_t j = 0; j < 4; j++)
        // {
        //   AK[cn * 16 + i * 4 + j] = 0;
        // }
      }
      else
      {
        d[I[i]] += AK[cn * 16 + i * 4 + i];
      }
    }
  }

  auto t1 = high_resolution_clock::now();
  for (size_t iterations = 0; iterations < max_iterations; iterations++)
  {
    // std::cout << "\rElement-by-Element Jacobi Iteration : " << iterations+1;

    std::memcpy(c, b, nV * sizeof(double));

    for (size_t cn = 0; cn < nC; cn++)
    {
      I[0] = mesh.Cells[cn].v0;
      I[1] = mesh.Cells[cn].v1;
      I[2] = mesh.Cells[cn].v2;
      I[3] = mesh.Cells[cn].v3;
      for (u_int8_t i = 0; i < 4; i++)
      {
        c[I[i]] -= AK[cn * 16 + i * 4 + (i + 1) % 4] * u[I[(i + 1) % 4]] +
                   AK[cn * 16 + i * 4 + (i + 2) % 4] * u[I[(i + 2) % 4]] +
                   AK[cn * 16 + i * 4 + (i + 3) % 4] * u[I[(i + 3) % 4]];
      }
    }

    // Update solution vector
    for (size_t i = 0; i < nV; i++)
    {
      u[i] = c[i] / d[i];
    }

    // TO DO: Check Convergance
  }
  auto t2 = high_resolution_clock::now();
  duration<double, std::milli> ms_double = t2 - t1;

  std::cout << "Jacobi finished after " << max_iterations << " iterations"
            << std::endl;
  std::cout << "Execution Time:" << ms_double.count() << "ms " << std::endl;

  // Calculating Error
  COO_array A_coo(nV, nV);
  assembleSparse(mesh, AK, isBoundary, &A_coo);

  CSR_array A(&A_coo);

  double errorA = calcError(&A, u, b);
  std::cout << "\n~~~Mean Absolute Error of  Jacobi Iterative Method: \n Err = "
            << errorA << std::endl;

  delete[] c, d, AK, u, b, isBoundary, vMarkers;
  return;
}

// Compute the number of Cells that each Vertex belongs
void getVertexDegree(uint *VertexDegree, const Mesh3D &mesh)
{
  for (size_t cn = 0; cn < mesh.Cells.size(); cn++)
  {
    VertexDegree[mesh.Cells[cn].v0]++;
    VertexDegree[mesh.Cells[cn].v1]++;
    VertexDegree[mesh.Cells[cn].v2]++;
    VertexDegree[mesh.Cells[cn].v3]++;
  }
}

void smoothLaplaceGaussSeidel(const Mesh3D &mesh,
                              const CityModel &citymodel,
                              const GridField2D &dtm,
                              const size_t max_iterations)
{
  info("Mesh Smoothing with Element-by-Element Gauss-Seidel Method");

  // const size_t max_iterations = MAX_ITERATIONS;
  const size_t nC = mesh.Cells.size();
  const size_t nV = mesh.Vertices.size();

  // Compute the number of Cells that each Vertex belongs
  uint *vertexDegree = new uint[nV]();
  uint *vDeg = new uint[nV]();
  getVertexDegree(vertexDegree, mesh);

  // Compute local Stiffness Matrices
  double *AK = new double[16 * nC];
  compute_transformation_matrix(AK, mesh);

  // Check Boundary Vertices
  bool *isBoundary = new bool[nV]();
  double *b = new double[nV]();
  checkBoundaryPoints(isBoundary, b, mesh);

  std::array<uint, 4> I = {0};

  // Diagonal d
  double *d = new double[nV]();

  // Initial guess for solution 0 TBD
  double *u = new double[nV]();
  std::memcpy(u, b, nV * sizeof(double));
  // Non-diagonal Elements sum
  double *c = new double[nV];

  auto t1 = high_resolution_clock::now();
  for (size_t iterations = 0; iterations < max_iterations; iterations++)
  {
    // std::cout << "\rElement-by-Element Jacobi Iteration : " << iterations+1;
    // memset(c, 0 ,nV * sizeof(double));
    std::memcpy(c, b, nV * sizeof(double));
    std::memcpy(vDeg, vertexDegree, nV * sizeof(uint));
    // std::copy(&b[0],&b[nV],c);

    for (size_t cn = 0; cn < nC; cn++)
    {
      I[0] = mesh.Cells[cn].v0;
      I[1] = mesh.Cells[cn].v1;
      I[2] = mesh.Cells[cn].v2;
      I[3] = mesh.Cells[cn].v3;
      for (u_int8_t i = 0; i < 4; i++)
      {
        c[I[i]] -= AK[cn * 16 + i * 4 + (i + 1) % 4] * u[I[(i + 1) % 4]] +
                   AK[cn * 16 + i * 4 + (i + 2) % 4] * u[I[(i + 2) % 4]] +
                   AK[cn * 16 + i * 4 + (i + 3) % 4] * u[I[(i + 3) % 4]];

        vDeg[I[i]]--;
        if (vDeg[I[i]] == 0)
        {
          u[I[i]] = c[I[i]] / d[I[i]];
        }
      }
    }

    // TO DO: Check Convergance
  }
  auto t2 = high_resolution_clock::now();
  duration<double, std::milli> ms_double = t2 - t1;
  std::cout << "GS finished after " << max_iterations << " iterations"
            << std::endl;
  std::cout << "Execution Time:" << ms_double.count() << "ms " << std::endl;

  COO_array A_coo(nV, nV);
  assembleSparse(mesh, AK, isBoundary, &A_coo);

  CSR_array A(&A_coo);
  double error = calcError(&A, u, b);
  std::cout << "\nMean Absolute Error of EbE GS Iterative Method: \n Err = "
            << error << std::endl;

  delete[] c, d, AK, u, b, isBoundary;
  delete[] vDeg;
  return;
}

#endif