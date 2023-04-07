#ifndef ASSEMBLED_H
#define ASSEMBLED_H

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

#include "error.hpp"
#include "sparse.hpp"
#include "stiffnessMatrix.hpp"

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

void assembleSparse(const Mesh3D &mesh,
                    double *A,
                    bool *isBoundary,
                    COO_array *assembled_A)
{
  size_t nC = mesh.Cells.size();
  size_t nV = mesh.Vertices.size();

  // double *assembled_A = new double[nV * nV];
  // std::cout << "\nNumber of Verticess: " << m.Vertices.size() << std::endl;
  std::array<uint, 4> I;

  for (size_t cn = 0; cn < nC; cn++)
  {
    I[0] = mesh.Cells[cn].v0;
    I[1] = mesh.Cells[cn].v1;
    I[2] = mesh.Cells[cn].v2;
    I[3] = mesh.Cells[cn].v3;

    for (size_t i = 0; i < 4; i++)
    {
      if (!isBoundary[I[i]])
      {
        assembled_A->add(I[i], I[0], A[16 * cn + i * 4 + 0]);
        assembled_A->add(I[i], I[1], A[16 * cn + i * 4 + 1]);
        assembled_A->add(I[i], I[2], A[16 * cn + i * 4 + 2]);
        assembled_A->add(I[i], I[3], A[16 * cn + i * 4 + 3]);
      }
    }
  }

  for (size_t i = 0; i < nV; i++)
  {
    if (isBoundary[i])
    {
      assembled_A->add(i, i, 1);
    }
  }
}

void assembled_GaussSeidel(const Mesh3D &mesh,
                           const CityModel &citymodel,
                           const GridField2D &dtm,
                           const size_t max_iterations)
{
  info("\nMesh Smoothing with Assembled Stiffness Matrix GS Method");

  const size_t nC = mesh.Cells.size();
  const size_t nV = mesh.Vertices.size();

  // Check Boundary Vertices
  bool *isBoundary = new bool[nV]();
  double *b = new double[nV]();

  checkBoundaryPoints(isBoundary, b, mesh);

  // Compute local Stiffness Matrices
  double *AK = new double[16 * nC];
  compute_transformation_matrix(AK, mesh);

  COO_array A_coo(nV, nV);
  assembleSparse(mesh, AK, isBoundary, &A_coo);

  CSR_array A(&A_coo);

  std::array<uint, 4> I = {0};
  // Diagonal d
  double *d = new double[nV]();

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
      }
      else
      {
        d[I[i]] += AK[cn * 16 + i * 4 + i];
      }
    }
  }

  // Initial guess for solution 0 TBD
  double *u_new = new double[nV]();
  double *u_old = new double[nV]();
  std::memcpy(u_new, b, nV * sizeof(double));

  // Non-diagonal Elements sum
  double c;

  auto t1 = high_resolution_clock::now();
  for (size_t iterations = 0; iterations < max_iterations; iterations++)
  {

    std::memcpy(u_old, u_new, nV * sizeof(double));
    c = 0;

    for (size_t i = 0; i < nV; i++)
    {
      c = 0;
      for (size_t j = A.rowPtr[i]; j < A.rowPtr[i + 1]; j++)
      {
        if (i != A.colIdx[j])
        {
          c += A.data[j] * u_new[A.colIdx[j]];
        }
      }
      u_new[i] = (b[i] - c) / d[i];
    }
  }
  auto t2 = high_resolution_clock::now();
  duration<double, std::milli> ms_double = t2 - t1;
  std::cout << "Execution Time:" << ms_double.count() << "ms " << std::endl;

  // std::memset(u_new,-10,nV * sizeof(double));

  double error = calcError(&A, u_new, b);
  std::cout << "\nMean Absolute Error of  GS Iterative Method: \n Err = "
            << error << std::endl;

  delete[] isBoundary, b, d, u_new, u_old, AK;
}

#endif
