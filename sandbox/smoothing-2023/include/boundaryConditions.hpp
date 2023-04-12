
#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

// #include "GridField.h"
// #include "LaplacianSmoother.h"
//#include "Mesh.h"
// #include "datamodel/CityModel.h"

using namespace DTCC;

void getVerticeMarkers(const Mesh3D &mesh, int *vMarkers)
{
  const size_t nC = mesh.Cells.size();
  const size_t nV = mesh.Vertices.size();

  for (size_t i = 0; i < nV; i++)
  {
    vMarkers[i] = -4;
  }

  std::array<uint, 4> I = {0};

  size_t k0 = 0;
  size_t k1 = 0;
  size_t k2 = 0;
  size_t k3 = 0;
  for (size_t cn = 0; cn < nC; cn++)
  {
    // Initializing Global Index for each cell
    I[0] = mesh.Cells[cn].v0;
    I[1] = mesh.Cells[cn].v1;
    I[2] = mesh.Cells[cn].v2;
    I[3] = mesh.Cells[cn].v3;

    const double z_mean = (mesh.Vertices[I[0]].z + mesh.Vertices[I[1]].z +
                           mesh.Vertices[I[2]].z + mesh.Vertices[I[3]].z) /
                          4;

    const int cellMarker = mesh.Markers[cn];

    if (cellMarker >= 0) // Building
    {
      for (size_t i = 0; i < 4; i++)
      {
        if (mesh.Vertices[I[i]].z > z_mean)
        {
          continue;
        }
        vMarkers[I[i]] = cellMarker;
      }
    }
    else if (cellMarker == -1) // Building Halo
    {
      for (size_t i = 0; i < 4; i++)
      {
        if (mesh.Vertices[I[i]].z > z_mean)
        {
          continue;
        }
        vMarkers[I[i]] = std::max(vMarkers[I[i]], -1);
      }
    }
    else if (cellMarker == -2) // Ground
    {
      for (size_t i = 0; i < 4; i++)
      {
        if (mesh.Vertices[I[i]].z > z_mean)
        {
          continue;
        }
        vMarkers[I[i]] = std::max(vMarkers[I[i]], -2);
      }
    }
    else if (cellMarker == -3) // Top
    {
      for (size_t i = 0; i < 4; i++)
      {
        if (mesh.Vertices[I[i]].z < z_mean)
        {
          continue;
        }
        vMarkers[I[i]] = std::max(vMarkers[I[i]], -3);
      }
    }
  }

  // Counting all vertice Markers
  for (size_t v = 0; v < nV; v++)
  {
    if (vMarkers[v] >= 0)
      k0++;
    else if (vMarkers[v] == -1)
      k1++;
    else if (vMarkers[v] == -2)
      k2++;
    else if (vMarkers[v] == -3)
      k3++;
  }

  int k4 = nV - (k0 + k1 + k2 + k3);
  std::cout << "Buildings k0 :" << k0 << std::endl
            << "Halos     k1 :" << k1 << std::endl
            << "Ground    k2 :" << k2 << std::endl
            << "Top       k3 :" << k3 << std::endl
            << "Neuman    k4 :" << k4 << std::endl
            << k0 + k1 + k2 + k3 + k4 << " / " << nV << std::endl;

  return;
}

void applyBC_load(const int *vMarkers,
                  const Mesh3D &mesh,
                  const CityModel &citymodel,
                  const GridField2D &elevation,
                  const double topHeight,
                  double *boundaryValues)
{
  const size_t nC = mesh.Cells.size();
  const size_t nV = mesh.Vertices.size();

  citymodel.BuildSearchTree();
  for (size_t i = 0; i < nV; i++)
  {
    const int verticeMarker = vMarkers[i];
    if (verticeMarker >= 0) // Building
    {
      boundaryValues[i] =
          citymodel.Buildings[verticeMarker].MaxHeight() - mesh.Vertices[i].z;
    }
    else if (verticeMarker == -1) // Building Halo
    {
      const Vector2D p(mesh.Vertices[i].x, mesh.Vertices[i].y);
      int buildingIndex = citymodel.FindBuilding(p);
      boundaryValues[i] =
          citymodel.Buildings[buildingIndex].MinHeight() - mesh.Vertices[i].z;
    }
    else if (verticeMarker == -2) // Ground
    {
      const Vector2D p(mesh.Vertices[i].x, mesh.Vertices[i].y);
      boundaryValues[i] = elevation(p) - mesh.Vertices[i].z;
    }
    else if (verticeMarker == -3) // Top
    {
      boundaryValues[i] = topHeight - mesh.Vertices[i].z;
    }
    else
    {
      boundaryValues[i] = 0;
    }
  }
}

void applyBC_stiffnessMat(const int *vMarkers,
                          const Mesh3D &mesh,
                          double *AK,
                          double *diag)
{
  const size_t nC = mesh.Cells.size();
  const size_t nV = mesh.Vertices.size();

  std::memset(diag, 0, nV * sizeof(double));
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
      if (vMarkers[I[i]] > -4)
      {
        diag[I[i]] = 1.0;
        // memset(&AK[cn * 16 + i * 4], 0, 4 * sizeof(double));
        for (uint8_t j = 0; j < 4; j++)
        {
          AK[cn * 16 + i * 4 + j] = 0;
        }
      }
      else
      {
        diag[I[i]] += AK[cn * 16 + i * 4 + i];
      }
    }
  }
}
#endif