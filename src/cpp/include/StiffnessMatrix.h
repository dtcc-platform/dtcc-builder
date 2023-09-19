// Copyright (C) 2023 George Spaias
// Licensed under the MIT License

#ifndef DTCC_STIFFNESS_MATRIX_H
#define DTCC_STIFFNESS_MATRIX_H

#include "model/Mesh.h"

namespace DTCC_BUILDER
{

/// Compute element matrix on given cell (tetrahedron)
//
// Input: vertices = array of length 12
// Represents flattened 4x3 matrix of vertex coordinates
//
// Output: A = array of length 16
// Represents flattened 4x4 matrix
void compute_element_matrix(double *A, const double *vertices)
{
  const double j_c4 = vertices[7] - vertices[1];
  const double j_c8 = vertices[11] - vertices[2];
  const double j_c5 = vertices[10] - vertices[1];
  const double j_c7 = vertices[8] - vertices[2];
  const double j_c0 = vertices[3] - vertices[0];
  const double j_c1 = vertices[6] - vertices[0];
  const double j_c6 = vertices[5] - vertices[2];
  const double j_c3 = vertices[4] - vertices[1];
  const double j_c2 = vertices[9] - vertices[0];

  alignas(32) double sp[80];
  sp[0] = j_c4 * j_c8;
  sp[1] = j_c5 * j_c7;
  sp[2] = sp[0] + -1 * sp[1];
  sp[3] = j_c0 * sp[2];
  sp[4] = j_c5 * j_c6;
  sp[5] = j_c3 * j_c8;
  sp[6] = sp[4] + -1 * sp[5];
  sp[7] = j_c1 * sp[6];
  sp[8] = sp[3] + sp[7];
  sp[9] = j_c3 * j_c7;
  sp[10] = j_c4 * j_c6;
  sp[11] = sp[9] + -1 * sp[10];
  sp[12] = j_c2 * sp[11];
  sp[13] = sp[8] + sp[12];
  sp[14] = sp[2] / sp[13];
  sp[15] = j_c3 * (-1 * j_c8);
  sp[16] = sp[4] + sp[15];
  sp[17] = sp[16] / sp[13];
  sp[18] = sp[11] / sp[13];
  sp[19] = sp[14] * sp[14];
  sp[20] = sp[14] * sp[17];
  sp[21] = sp[18] * sp[14];
  sp[22] = sp[17] * sp[17];
  sp[23] = sp[18] * sp[17];
  sp[24] = sp[18] * sp[18];
  sp[25] = j_c2 * j_c7;
  sp[26] = j_c8 * (-1 * j_c1);
  sp[27] = sp[25] + sp[26];
  sp[28] = sp[27] / sp[13];
  sp[29] = j_c0 * j_c8;
  sp[30] = j_c6 * (-1 * j_c2);
  sp[31] = sp[29] + sp[30];
  sp[32] = sp[31] / sp[13];
  sp[33] = j_c1 * j_c6;
  sp[34] = j_c0 * j_c7;
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
  sp[49] = j_c1 * j_c5;
  sp[50] = j_c2 * j_c4;
  sp[51] = sp[49] + -1 * sp[50];
  sp[52] = sp[51] / sp[13];
  sp[53] = j_c2 * j_c3;
  sp[54] = j_c0 * j_c5;
  sp[55] = sp[53] + -1 * sp[54];
  sp[56] = sp[55] / sp[13];
  sp[57] = j_c0 * j_c4;
  sp[58] = j_c1 * j_c3;
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

// Stiffness matrix (unassembled)
class StiffnessMatrix
{
public:
  // Data as a contiguous array
  double *_data;

  // Matrix shape (number of cells, 4, 4)
  std::array<size_t, 3> shape;

  // Diagonal elements
  std::vector<double> diagonal;

  // Constructor
  StiffnessMatrix(const VolumeMesh &volume_mesh)
      : _data(0), diagonal(volume_mesh.vertices.size(), 0)
  {
    // Set shape
    shape[0] = volume_mesh.cells.size();
    shape[1] = 4;
    shape[2] = 4;

    // Compute stiffness matrix
    _data = new double[shape[0] * shape[1] * shape[2]];
    compute_stiffness_matrix(volume_mesh);

    // Fill diagonal vector with appropriate values
    for (size_t c = 0; c < volume_mesh.cells.size(); c++)
    {
      diagonal[volume_mesh.cells[c].v0] += _data[c * 16 + 0 * 4 + 0];
      diagonal[volume_mesh.cells[c].v1] += _data[c * 16 + 1 * 4 + 1];
      diagonal[volume_mesh.cells[c].v2] += _data[c * 16 + 2 * 4 + 2];
      diagonal[volume_mesh.cells[c].v3] += _data[c * 16 + 3 * 4 + 3];
    }
  }

  // Destructor
  ~StiffnessMatrix() { delete[] _data; }

  // Access matrix elements
  inline double &operator()(unsigned cell, unsigned row, unsigned col) const
  {
    return _data[16 * cell + 4 * row + col];
  }

private:
  // Compute stiffness matrix
  void compute_stiffness_matrix(const VolumeMesh &volume_mesh)
  {
    // Array for vertex coordinates
    double vertices[12];

    // Iterate over cells
    for (size_t c = 0; c < volume_mesh.cells.size(); c++)
    {
      // Set vertex coordinates
      vertices[0] = volume_mesh.vertices[volume_mesh.cells[c].v0].x;
      vertices[1] = volume_mesh.vertices[volume_mesh.cells[c].v0].y;
      vertices[2] = volume_mesh.vertices[volume_mesh.cells[c].v0].z;
      vertices[3] = volume_mesh.vertices[volume_mesh.cells[c].v1].x;
      vertices[4] = volume_mesh.vertices[volume_mesh.cells[c].v1].y;
      vertices[5] = volume_mesh.vertices[volume_mesh.cells[c].v1].z;
      vertices[6] = volume_mesh.vertices[volume_mesh.cells[c].v2].x;
      vertices[7] = volume_mesh.vertices[volume_mesh.cells[c].v2].y;
      vertices[8] = volume_mesh.vertices[volume_mesh.cells[c].v2].z;
      vertices[9] = volume_mesh.vertices[volume_mesh.cells[c].v3].x;
      vertices[10] = volume_mesh.vertices[volume_mesh.cells[c].v3].y;
      vertices[11] = volume_mesh.vertices[volume_mesh.cells[c].v3].z;

      // Compute element matrix
      compute_element_matrix(_data + c * 16, vertices);
    }
  }
};

} // namespace DTCC_BUILDER

#endif
