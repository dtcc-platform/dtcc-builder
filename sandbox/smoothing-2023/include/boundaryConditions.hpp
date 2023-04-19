
#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "GridField.h"
// #include "LaplacianSmoother.h"
// #include "Mesh.h"
// #include "datamodel/CityModel.h"

#include "stiffnessMatrix.hpp"

using namespace DTCC;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Boundary Conditions Class (WIP)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class BoundaryConditions
{
private:
  Mesh3D &_mesh;

  const CityModel &_citymodel;

  const GridField2D &_dtm;

  const double topHeight;

  const bool fixBuildings;

public:
  std::vector<int> vMarkers;

  std::vector<double> values;

  BoundaryConditions(Mesh3D &mesh,
                     const CityModel &cityModel,
                     const GridField2D &dtm,
                     const double topHeight,
                     const bool fixBuildings);

  ~BoundaryConditions();

  void computeVerticeMarkers();

  void computeBoundaryValues();

  void apply(std::vector<double> &b);

  void apply(stiffnessMatrix &A);
};

BoundaryConditions::BoundaryConditions(Mesh3D &mesh,
                                       const CityModel &cityModel,
                                       const GridField2D &dtm,
                                       const double topHeight,
                                       const bool fixBuildings)
    : _mesh(mesh), _citymodel(cityModel), _dtm(dtm), topHeight(topHeight),
      vMarkers(mesh.Vertices.size(), 0), values(mesh.Vertices.size(), 0.0),
      fixBuildings(fixBuildings)
{
  const size_t nC = _mesh.Cells.size();
  const size_t nV = _mesh.Vertices.size();

  info("Translating Cell Markers to Vertice Markers");
  computeVerticeMarkers();

  info("Computing Boundary Values");
  computeBoundaryValues();
}

BoundaryConditions::~BoundaryConditions() {}

void BoundaryConditions::computeVerticeMarkers()
{
  const size_t nC = _mesh.Cells.size();
  const size_t nV = _mesh.Vertices.size();

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
    I[0] = _mesh.Cells[cn].v0;
    I[1] = _mesh.Cells[cn].v1;
    I[2] = _mesh.Cells[cn].v2;
    I[3] = _mesh.Cells[cn].v3;

    const double z_mean = (_mesh.Vertices[I[0]].z + _mesh.Vertices[I[1]].z +
                           _mesh.Vertices[I[2]].z + _mesh.Vertices[I[3]].z) /
                          4;

    const int cellMarker = _mesh.Markers[cn];

    if (cellMarker >= 0 && fixBuildings) // Building
    {
      for (size_t i = 0; i < 4; i++)
      {
        if (_mesh.Vertices[I[i]].z > z_mean)
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
        if (_mesh.Vertices[I[i]].z > z_mean)
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
        if (_mesh.Vertices[I[i]].z > z_mean)
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
        if (_mesh.Vertices[I[i]].z < z_mean)
        {
          continue;
        }
        vMarkers[I[i]] = std::max(vMarkers[I[i]], -3);
      }
    }
  }

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
  std::cout << "Buildings      k0 :" << k0 << std::endl
            << "Building Halos k1 :" << k1 << std::endl
            << "Ground         k2 :" << k2 << std::endl
            << "Top            k3 :" << k3 << std::endl
            << "Neumann        k4 :" << k4 << std::endl
            << k0 + k1 + k2 + k3 + k4 << "  " << nV << std::endl;

  return;
}

void BoundaryConditions::computeBoundaryValues()
{
  const size_t nC = _mesh.Cells.size();
  const size_t nV = _mesh.Vertices.size();

  // TODO: Check if Search tree has already been built
  //_citymodel.BuildSearchTree(true);

  for (size_t i = 0; i < nV; i++)
  {
    const int verticeMarker = vMarkers[i];
    if (verticeMarker >= 0 && fixBuildings) // Building
    {
      values[i] =
          _citymodel.Buildings[verticeMarker].MaxHeight() - _mesh.Vertices[i].z;
    }
    else if (verticeMarker == -1) // Building Halo
    {
      const Vector2D p(_mesh.Vertices[i].x, _mesh.Vertices[i].y);
      int buildingIndex = _citymodel.FindBuilding(p);
      values[i] =
          _citymodel.Buildings[buildingIndex].MinHeight() - _mesh.Vertices[i].z;
    }
    else if (verticeMarker == -2) // Ground
    {
      const Vector2D p(_mesh.Vertices[i].x, _mesh.Vertices[i].y);
      values[i] = _dtm(p) - _mesh.Vertices[i].z;
    }
    else if (verticeMarker == -3) // Top
    {
      values[i] = topHeight - _mesh.Vertices[i].z;
    }
    else
    {
      values[i] = 0;
    }
  }
}

// Apply Boundary Conditions on Load vector
void BoundaryConditions::apply(std::vector<double> &b)
{
  info("Applying Boundary Conditions on Load vector b");
  b = values;
}

// Apply Boundary Conditions on Stiffness Matrix
void BoundaryConditions::apply(stiffnessMatrix &A)
{
  info("Applying Boundary Conditions on Stiffness Matrix AK");

  std::array<uint, 4> I = {0};

  for (size_t cn = 0; cn < A.shape[0]; cn++)
  {
    // Global Index for each cell
    I[0] = _mesh.Cells[cn].v0;
    I[1] = _mesh.Cells[cn].v1;
    I[2] = _mesh.Cells[cn].v2;
    I[3] = _mesh.Cells[cn].v3;

    for (size_t i = 0; i < A.shape[1]; i++)
    {
      if (vMarkers[I[i]] > -4)
      {
        A.diagonal[I[i]] = 1.0;
        for (size_t j = 0; j < A.shape[2]; j++)
        {
          A(cn, i, j) = 0;
        }
      }
    }
  }
}

#endif