
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
  // Vertex Boundary Markers:
  // -4 : Neumann Vertices
  // -3 : Top Boundary Vertices
  // -2 : Ground Boundary Vertices
  // -1 : Building Halos Boundary Vertices
  std::vector<int> vMarkers;

  // Boundary Values
  std::vector<double> values;

  // Building Polygon Centroids
  std::vector<Vector2D> BuildingCentroids;

  // Elevation for halo vertices based on min Cell elevation
  std::vector<double> HaloElevations;

  BoundaryConditions(Mesh3D &mesh,
                     const CityModel &cityModel,
                     const GridField2D &dtm,
                     const double topHeight,
                     const bool fixBuildings);

  ~BoundaryConditions();

  void computeVerticeMarkers();

  void computeBoundaryValues();

  void computeBuildingCentroids();

  int FindAdjacentBuilding(const Vector2D &p);

  void HaloExpressionDEM();

  void apply(std::vector<double> &b);

  void apply(stiffnessMatrix &A);
};

BoundaryConditions::BoundaryConditions(Mesh3D &mesh,
                                       const CityModel &cityModel,
                                       const GridField2D &dtm,
                                       const double topHeight,
                                       const bool fixBuildings)
    : _mesh(mesh), _citymodel(cityModel), _dtm(dtm), topHeight(topHeight),
      vMarkers(mesh.Vertices.size(), -4), values(mesh.Vertices.size(), 0.0),
      fixBuildings(fixBuildings),
      HaloElevations(mesh.Vertices.size(), std::numeric_limits<double>::max())
{
  const size_t nC = _mesh.Cells.size();
  const size_t nV = _mesh.Vertices.size();

  info("Translating Cell Markers to Vertice Markers");
  computeVerticeMarkers();

  info("Computing Boundary Values");
  computeBoundaryValues();
}

BoundaryConditions::~BoundaryConditions() {}

// Compute Vertex Boundary Markers based on Cell Boundary Markers
void BoundaryConditions::computeVerticeMarkers()
{
  const size_t nC = _mesh.Cells.size();
  const size_t nV = _mesh.Vertices.size();
  std::array<uint, 4> I = {0};

  size_t k0 = 0;
  size_t k1 = 0;
  size_t k2 = 0;
  size_t k3 = 0;

  const double domainMin = _dtm.Min();
  const double domainMean = _dtm.Mean();
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
    const double BuildingMaxHeight =
        _citymodel.Buildings[cellMarker].MaxHeight();
    const double BuildingMinHeight =
        _citymodel.Buildings[cellMarker].MinHeight();
    if (cellMarker >= 0 && fixBuildings) // Building
    {
      for (size_t i = 0; i < 4; i++)
      {
        // Test: This if can be removed completely
        /*   if (_mesh.Vertices[I[i]].z > (BuildingMaxHeight))
          {
            std::cout << I[i] << ") Vertex z: " << _mesh.Vertices[I[i]].z
                      << " B_ID: " << cellMarker
                      << " Building Height: " << BuildingMaxHeight << std::endl;
            // continue;
          } */
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

// Compute Values for Boundary Vertices
void BoundaryConditions::computeBoundaryValues()
{
  const size_t nC = _mesh.Cells.size();
  const size_t nV = _mesh.Vertices.size();

  // TODO: Check if Search tree has already been built
  //_citymodel.BuildSearchTree(true,0.0);

  // Min adjacent Building Height Expression for Halos
  computeBuildingCentroids();

  // DEM expression for Halos
  HaloExpressionDEM();

  std::size_t k1 = 0;
  for (size_t i = 0; i < nV; i++)
  {
    const int verticeMarker = vMarkers[i];
    if (verticeMarker >= 0) //  && fixBuildings Building
    {
      values[i] =
          _citymodel.Buildings[verticeMarker].MaxHeight() - _mesh.Vertices[i].z;
    }
    else if (verticeMarker == -1) // Building Halo
    {

      // Halo Min Building Height Expression
      const Vector2D p(_mesh.Vertices[i].x, _mesh.Vertices[i].y);
      const int buildingIndex = FindAdjacentBuilding(p);

      if (buildingIndex == -1)
      {
        // std::cout << "Vertex " << i << " is not adjacent to any building...
        // "<< std::endl;
        values[i] = _dtm(p) - _mesh.Vertices[i].z;
        k1++;
      }
      else
      {
        values[i] = _citymodel.Buildings[buildingIndex].MinHeight() -
                    _mesh.Vertices[i].z;
      }

      // Halo DEM Expression
      // values[i] = HaloElevations[i] - _mesh.Vertices[i].z;
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

  std::cout << "Didnt find Building Index for " << k1 << " Halo Vertices"
            << std::endl;
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

// Finds Building Polygon Centroids
void BoundaryConditions::computeBuildingCentroids()
{
  const std::size_t nV = _mesh.Vertices.size();

  BuildingCentroids.resize(_citymodel.Buildings.size());

  for (size_t i = 0; i < _citymodel.Buildings.size(); i++)
  {
    Vector2D p(0, 0);
    Polygon fp = _citymodel.Buildings[i].Footprint;

    for (auto vertex : fp.Vertices)
    {
      p += Vector2D(vertex);
    }
    p = p / static_cast<double>(fp.Vertices.size());
    BuildingCentroids[i] = p;
  }
}

// Returns the index of the closest Building that is closest to the input point
int BoundaryConditions::FindAdjacentBuilding(const Vector2D &p)
{
  int building_idx = -1;
  double minDistance = MAXFLOAT;
  const size_t numOfBuildings = _citymodel.Buildings.size();

  for (size_t i = 0; i < numOfBuildings; i++)
  {
    Vector2D q = p - BuildingCentroids[i];
    double distance = q.SquaredMagnitude();

    if (distance < minDistance)
    {
      minDistance = distance;
      building_idx = i;
    }
  }

  return building_idx;
}

// Compute Halo Boundary vertices elevation based on
// the min elevation in the cell containing them
void BoundaryConditions::HaloExpressionDEM()
{
  const std::size_t nC = _mesh.Cells.size();
  const std::size_t nV = _mesh.Vertices.size();

  std::array<uint, 4> I = {0};

  for (size_t cn = 0; cn < nC; cn++)
  {
    I[0] = _mesh.Cells[cn].v0;
    I[1] = _mesh.Cells[cn].v1;
    I[2] = _mesh.Cells[cn].v2;
    I[3] = _mesh.Cells[cn].v3;

    double z_min = std::numeric_limits<double>::max();

    for (size_t i = 0; i < 4; i++)
    {
      const Vector2D p(_mesh.Vertices[I[i]].x, _mesh.Vertices[I[i]].y);
      const double z = _dtm(p);

      z_min = std::min(z_min, z);
    }

    HaloElevations[I[0]] = std::min(HaloElevations[I[0]], z_min);
    HaloElevations[I[1]] = std::min(HaloElevations[I[1]], z_min);
    HaloElevations[I[2]] = std::min(HaloElevations[I[2]], z_min);
    HaloElevations[I[3]] = std::min(HaloElevations[I[3]], z_min);
  }
}

#endif