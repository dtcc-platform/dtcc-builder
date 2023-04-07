#ifndef ERROR_HPP
#define ERROR_HPP

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "Mesh.h"
#include "sparse.hpp"

using namespace DTCC;

// Returns Error err= b-A*u for our linear System solution
double calcError(double *A, double *u, double *b, size_t nV)
{
  double *error = new double[nV]();

  if (nV == 0)
    return -1;
  for (size_t i = 0; i < nV; i++)
  {
    for (size_t j = 0; j < nV; j++)
    {
      error[i] += A[i * nV + j] * u[j];
    }
    error[i] -= b[i];
  }

  long double sumOfSquares = 0;
  for (size_t i = 0; i < nV; i++)
  {
    sumOfSquares += error[i] * error[i];
  }

  delete[] error;
  return std::sqrt(sumOfSquares);
}

double calcError(
    const Mesh3D *mesh, double *AK, double *u, double *b, bool *isBoundary)
{
  const size_t nC = mesh->Cells.size();
  const size_t nV = mesh->Vertices.size();

  double *error = new double[nV]();
  std::memset(error, 0, nV * sizeof(double));
  std::array<uint, 4> I;

  for (size_t cn = 0; cn < nC; cn++)
  {
    I[0] = mesh->Cells[cn].v0;
    I[1] = mesh->Cells[cn].v1;
    I[2] = mesh->Cells[cn].v2;
    I[3] = mesh->Cells[cn].v3;

    for (uint i = 0; i < 4; i++)
    {
      for (uint8_t j = 0; j < 4; j++)
      {
        error[I[i]] += AK[cn * 16 + i * 4 + j] * u[I[j]];
      }
    }
  }

  double sumOfSquares = 0;
  for (size_t i = 0; i < nV; i++)
  {
    if (isBoundary[i])
    {
      error[i] = u[i] - b[i];
    }
    else
    {
      error[i] -= b[i];
    }
    sumOfSquares += error[i] * error[i];
  }

  delete[] error;
  return std::sqrt(sumOfSquares);
}

double calcError(CSR_array *A, double *u, double *b)
{
  const size_t nV = A->shape[0];
  if (nV == 0)
  {
    return -1.0; // or throw an exception
  }

  double sumOfSquares = 0.0;
  for (size_t i = 0; i < nV; i++)
  {
    double error = -b[i];
    const size_t rowStart = A->rowPtr[i];
    const size_t rowEnd = A->rowPtr[i + 1];
    for (size_t k = rowStart; k < rowEnd; k++)
    {
      error += A->data[k] * u[A->colIdx[k]];
    }
    sumOfSquares += error * error;
  }
  return sqrt(sumOfSquares);
}

// {
//   const size_t nV = A->shape[0];
//   std::vector<long double> error(A->shape[0], 0.0);

//   if (nV == 0)
//     return -1;
//   for (size_t i = 0; i < nV; i++)
//   {
//     error[i] = 0;
//     for (size_t k = A->rowPtr[i]; k < A->rowPtr[i + 1]; k++)
//     {
//       error[i] += A->data[k] * u[A->colIdx[k]];
//     }
//     error[i] -= b[i];
//   }

//   long double sumOfSquares = 0;
//   for (size_t i = 0; i < nV; i++)
//   {
//     sumOfSquares += error[i] * error[i];
//   }
//   return std::sqrt(sumOfSquares);
// }

#endif