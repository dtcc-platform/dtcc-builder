// Copyright (C) 2020-2021 Anders Logg, Anton J Olsson
// Licensed under the MIT License

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Needs to come before JSON (nlohmann) include because of sloppy
// namespacing in VTK (typedef detail)...

// DTCC includes
#include "JSON.h"
#include "Mesh.h"
#include "Timer.h"

using namespace DTCC_BUILDER;

void Help() { error("Usage: dtcc-sandbox CityMesh.json"); }

void compute_element_matrix(double *A, const double *vertices);

void compute_transformation_matrix(double *A, Mesh3D *m);

int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (argc != 2)
  {
    Help();
    return 1;
  }

  // Read parameters
  Parameters p;
  Mesh3D m;
  JSON::Read(m, argv[1]);
  info(m);

  double *AA = new double[16 * m.Cells.size()];
  std::cout << "\nComputing Element matrix for all Cells\n";

  Timer t = Timer("Transformation Matrix Computation");
  compute_transformation_matrix(AA, &m);
  t.Stop();
  std::cout << "\nElapsed time for A matrix Computation:" << std::endl;
  t.Print();

  // std::ofstream outfile ("test.txt");
  // for (size_t i = 0; i < 16 * m.Cells.size(); i++)
  // {
  //   outfile << AA[i] << std::endl;
  // }
  // outfile.close();

  std::cout << "\nNumber of Cells: " << m.Cells.size() << std::endl;
  delete AA;
  return 0;
}

// Output: A = array of length 16
// Represents flattened 4x4 matrix
//
// Input: vertices = array of length 12
// Represents flattened 4x3 matrix of vertex coordinates

void compute_element_matrix(double *A, const double *vertices)
{
  const double J_c4 = vertices[7] - vertices[1];
  const double J_c8 = vertices[11] - vertices[2];
  const double J_c5 = vertices[10] - vertices[1];
  const double J_c7 = vertices[8] - vertices[2];
  const double J_c0 = vertices[3] - vertices[0];
  const double J_c1 = vertices[6] - vertices[0];
  const double J_c6 = vertices[5] - vertices[2];
  const double J_c3 = vertices[4] - vertices[1];
  const double J_c2 = vertices[9] - vertices[0];

  alignas(32) double sp[80];
  sp[0] = J_c4 * J_c8;
  sp[1] = J_c5 * J_c7;
  sp[2] = sp[0] + -1 * sp[1];
  sp[3] = J_c0 * sp[2];
  sp[4] = J_c5 * J_c6;
  sp[5] = J_c3 * J_c8;
  sp[6] = sp[4] + -1 * sp[5];
  sp[7] = J_c1 * sp[6];
  sp[8] = sp[3] + sp[7];
  sp[9] = J_c3 * J_c7;
  sp[10] = J_c4 * J_c6;
  sp[11] = sp[9] + -1 * sp[10];
  sp[12] = J_c2 * sp[11];
  sp[13] = sp[8] + sp[12];
  sp[14] = sp[2] / sp[13];
  sp[15] = J_c3 * (-1 * J_c8);
  sp[16] = sp[4] + sp[15];
  sp[17] = sp[16] / sp[13];
  sp[18] = sp[11] / sp[13];
  sp[19] = sp[14] * sp[14];
  sp[20] = sp[14] * sp[17];
  sp[21] = sp[18] * sp[14];
  sp[22] = sp[17] * sp[17];
  sp[23] = sp[18] * sp[17];
  sp[24] = sp[18] * sp[18];
  sp[25] = J_c2 * J_c7;
  sp[26] = J_c8 * (-1 * J_c1);
  sp[27] = sp[25] + sp[26];
  sp[28] = sp[27] / sp[13];
  sp[29] = J_c0 * J_c8;
  sp[30] = J_c6 * (-1 * J_c2);
  sp[31] = sp[29] + sp[30];
  sp[32] = sp[31] / sp[13];
  sp[33] = J_c1 * J_c6;
  sp[34] = J_c0 * J_c7;
  sp[35] = sp[33] + -1 * sp[34];
  sp[36] = sp[35] / sp[13];
  sp[37] = sp[28] * sp[28];
  sp[38] = sp[28] * sp[32];
  sp[39] = sp[28] * sp[36];
  sp[40] = sp[32] * sp[32];
  sp[41] = sp[32] * sp[36];
  sp[42] = sp[36] * sp[36];
  sp[43] = sp[37] + sp[19];
  sp[44] = sp[38] + sp[20];
  sp[45] = sp[39] + sp[21];
  sp[46] = sp[40] + sp[22];
  sp[47] = sp[41] + sp[23];
  sp[48] = sp[24] + sp[42];
  sp[49] = J_c1 * J_c5;
  sp[50] = J_c2 * J_c4;
  sp[51] = sp[49] + -1 * sp[50];
  sp[52] = sp[51] / sp[13];
  sp[53] = J_c2 * J_c3;
  sp[54] = J_c0 * J_c5;
  sp[55] = sp[53] + -1 * sp[54];
  sp[56] = sp[55] / sp[13];
  sp[57] = J_c0 * J_c4;
  sp[58] = J_c1 * J_c3;
  sp[59] = sp[57] + -1 * sp[58];
  sp[60] = sp[59] / sp[13];
  sp[61] = sp[52] * sp[52];
  sp[62] = sp[52] * sp[56];
  sp[63] = sp[60] * sp[52];
  sp[64] = sp[56] * sp[56];
  sp[65] = sp[60] * sp[56];
  sp[66] = sp[60] * sp[60];
  sp[67] = sp[43] + sp[61];
  sp[68] = sp[44] + sp[62];
  sp[69] = sp[45] + sp[63];
  sp[70] = sp[46] + sp[64];
  sp[71] = sp[47] + sp[65];
  sp[72] = sp[48] + sp[66];
  sp[73] = std::abs(sp[13]);
  sp[74] = sp[67] * sp[73];
  sp[75] = sp[68] * sp[73];
  sp[76] = sp[69] * sp[73];
  sp[77] = sp[70] * sp[73];
  sp[78] = sp[71] * sp[73];
  sp[79] = sp[72] * sp[73];

  A[0] = 0.1666666666666667 * sp[74] + 0.1666666666666667 * sp[75] +
         0.1666666666666667 * sp[76] + 0.1666666666666667 * sp[75] +
         0.1666666666666667 * sp[77] + 0.1666666666666667 * sp[78] +
         0.1666666666666667 * sp[76] + 0.1666666666666667 * sp[78] +
         0.1666666666666667 * sp[79];
  A[1] = -0.1666666666666667 * sp[74] + -0.1666666666666667 * sp[75] +
         -0.1666666666666667 * sp[76];
  A[2] = -0.1666666666666667 * sp[75] + -0.1666666666666667 * sp[77] +
         -0.1666666666666667 * sp[78];
  A[3] = -0.1666666666666667 * sp[76] + -0.1666666666666667 * sp[78] +
         -0.1666666666666667 * sp[79];
  A[4] = -0.1666666666666667 * sp[74] + -0.1666666666666667 * sp[75] +
         -0.1666666666666667 * sp[76];
  A[5] = 0.1666666666666667 * sp[74];
  A[6] = 0.1666666666666667 * sp[75];
  A[7] = 0.1666666666666667 * sp[76];
  A[8] = -0.1666666666666667 * sp[75] + -0.1666666666666667 * sp[77] +
         -0.1666666666666667 * sp[78];
  A[9] = 0.1666666666666667 * sp[75];
  A[10] = 0.1666666666666667 * sp[77];
  A[11] = 0.1666666666666667 * sp[78];
  A[12] = -0.1666666666666667 * sp[76] + -0.1666666666666667 * sp[78] +
          -0.1666666666666667 * sp[79];
  A[13] = 0.1666666666666667 * sp[76];
  A[14] = 0.1666666666666667 * sp[78];
  A[15] = 0.1666666666666667 * sp[79];
}

/*
 * Output: A = array of length 16 x Num of Cells in 3D Mesh
 * Represents flattened 4x4xNum of Cells matrix
 *
 * Input: vertices = Mesh3D Mesh
 * Represents flattened 4x3 matrix of vertex coordinates
 */

void compute_transformation_matrix(double *A, Mesh3D *m)
{
  const uint8_t size_A = 16;
  double *V = new double[12];
  size_t n = m->Cells.size();
  size_t i;
  for (i = 0; i < n; i++)
  {

    V[0] = m->Vertices[m->Cells[i].v0].x;
    V[1] = m->Vertices[m->Cells[i].v0].y;
    V[2] = m->Vertices[m->Cells[i].v0].z;

    V[3] = m->Vertices[m->Cells[i].v1].x;
    V[4] = m->Vertices[m->Cells[i].v1].y;
    V[5] = m->Vertices[m->Cells[i].v1].z;

    V[6] = m->Vertices[m->Cells[i].v2].x;
    V[7] = m->Vertices[m->Cells[i].v2].y;
    V[8] = m->Vertices[m->Cells[i].v2].z;

    V[9] = m->Vertices[m->Cells[i].v3].x;
    V[10] = m->Vertices[m->Cells[i].v3].y;
    V[11] = m->Vertices[m->Cells[i].v3].z;

    compute_element_matrix(A + size_A * i, V);
  }
  std::cout << i << std::endl;
  delete V;
}