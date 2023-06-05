// Copyright (C) 2023 George Spaias
// Licensed under the MIT License

#ifndef DTCC_BOUNDARY_CONDITIONS_H
#define DTCC_BOUNDARY_CONDITIONS_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "GridField.h"
#include "StiffnessMatrix.h"

namespace DTCC_BUILDER
{

class BoundaryConditions
{
public:
  // Vertex Boundary Markers:
  // -4 : Neumann Vertices
  // -3 : Top Boundary Vertices
  // -2 : Ground Boundary Vertices
  // -1 : Building Halos Boundary Vertices
  std::vector<int> vertex_markers;

  // Boundary Values
  std::vector<double> values;

  // Building Polygon Centroids
  std::vector<Vector2D> building_centroids;

  // Elevation for halo vertices based on min Cell elevation
  std::vector<double> halo_elevations;

  // Constructor
  BoundaryConditions(Mesh3D &mesh,
                     const CityModel &cityModel,
                     const GridField2D &dtm,
                     const double top_height,
                     const bool fix_buildings)
      : _mesh(mesh), _citymodel(cityModel), _dtm(dtm), top_height(top_height),
        vertex_markers(mesh.Vertices.size(), -4),
        values(mesh.Vertices.size(), 0.0), fix_buildings(fix_buildings),
        halo_elevations(mesh.Vertices.size(),
                        std::numeric_limits<double>::max())
  {
    // Compute vertex markers
    compute_vertex_markers();

    // Compute boundary values
    compute_boundary_values();
  }

  // Destructor
  ~BoundaryConditions() {}

  // Compute Vertex Boundary Markers based on Cell Boundary Markers
  void compute_vertex_markers()
  {
    info("BoundaryConditions: Computing vertex markers from cell markers");

    std::array<uint, 4> I = {0};

    size_t k0 = 0;
    size_t k1 = 0;
    size_t k2 = 0;
    size_t k3 = 0;

    for (size_t c = 0; c < _mesh.Cells.size(); c++)
    {
      // Initializing Global Index for each cell
      I[0] = _mesh.Cells[c].v0;
      I[1] = _mesh.Cells[c].v1;
      I[2] = _mesh.Cells[c].v2;
      I[3] = _mesh.Cells[c].v3;

      const double z_mean = (_mesh.Vertices[I[0]].z + _mesh.Vertices[I[1]].z +
                             _mesh.Vertices[I[2]].z + _mesh.Vertices[I[3]].z) /
                            4;

      const int cell_marker = _mesh.Markers[c];
      const double BuildingMaxHeight =
          _citymodel.Buildings[cell_marker].MaxHeight();
      const double BuildingMinHeight =
          _citymodel.Buildings[cell_marker].MinHeight();
      if (cell_marker >= 0 && fix_buildings) // Building
      {
        for (size_t i = 0; i < 4; i++)
        {
          if (_mesh.Vertices[I[i]].z > z_mean)
          {
            continue;
          }
          vertex_markers[I[i]] = cell_marker;
        }
      }
      else if (cell_marker == -1) // Building Halo
      {
        for (size_t i = 0; i < 4; i++)
        {
          if (_mesh.Vertices[I[i]].z > z_mean)
          {
            continue;
          }
          vertex_markers[I[i]] = std::max(vertex_markers[I[i]], -1);
        }
      }
      else if (cell_marker == -2) // Ground
      {
        for (size_t i = 0; i < 4; i++)
        {
          if (_mesh.Vertices[I[i]].z > z_mean)
          {
            continue;
          }
          vertex_markers[I[i]] = std::max(vertex_markers[I[i]], -2);
        }
      }
      else if (cell_marker == -3) // Top
      {
        for (size_t i = 0; i < 4; i++)
        {
          if (_mesh.Vertices[I[i]].z < z_mean)
          {
            continue;
          }
          vertex_markers[I[i]] = std::max(vertex_markers[I[i]], -3);
        }
      }
    }

    for (size_t v = 0; v < _mesh.Vertices.size(); v++)
    {
      if (vertex_markers[v] >= 0)
        k0++;
      else if (vertex_markers[v] == -1)
        k1++;
      else if (vertex_markers[v] == -2)
        k2++;
      else if (vertex_markers[v] == -3)
        k3++;
    }

    return;
  }

  // Compute boundary values
  void compute_boundary_values()
  {
    info("BoundaryConditions: Computing boundary values");

    // TODO: Check if Search tree has already been built
    //_citymodel.BuildSearchTree(true,0.0);

    // Compute building centroids
    compute_building_centroids();

    // Compute halo elevations
    compute_halo_elevations();

    for (size_t i = 0; i < _mesh.Vertices.size(); i++)
    {
      const int vertex_marker = vertex_markers[i];
      if (vertex_marker >= 0) //  && fix_buildings Building
      {
        values[i] = _citymodel.Buildings[vertex_marker].MaxHeight() -
                    _mesh.Vertices[i].z;
      }
      else if (vertex_marker == -1) // Building Halo
      {
        values[i] = halo_elevations[i] - _mesh.Vertices[i].z;
      }
      else if (vertex_marker == -2) // Ground
      {
        const Vector2D p(_mesh.Vertices[i].x, _mesh.Vertices[i].y);
        values[i] = _dtm(p) - _mesh.Vertices[i].z;
      }
      else if (vertex_marker == -3) // Top
      {
        values[i] = top_height - _mesh.Vertices[i].z;
      }
      else
      {
        values[i] = 0;
      }
    }
  }

  // Apply boundary conditions to stiffness matrix
  void apply(StiffnessMatrix &A)
  {
    info(
        "BoundaryConditions: Applying boundary conditions to stiffness matrix");

    std::array<uint, 4> I = {0};

    for (size_t c = 0; c < A.shape[0]; c++)
    {
      // Global Index for each cell
      I[0] = _mesh.Cells[c].v0;
      I[1] = _mesh.Cells[c].v1;
      I[2] = _mesh.Cells[c].v2;
      I[3] = _mesh.Cells[c].v3;

      for (size_t i = 0; i < A.shape[1]; i++)
      {
        if (vertex_markers[I[i]] > -4)
        {
          A.diagonal[I[i]] = 1.0;
          for (size_t j = 0; j < A.shape[2]; j++)
          {
            A(c, i, j) = 0;
          }
        }
      }
    }
  }

  // Apply boundary conditions to load vector
  void apply(std::vector<double> &b)
  {
    info("BoundaryConditions: Applying boundary conditions to load vector");
    b = values;
  }

private:
  Mesh3D &_mesh;

  const CityModel &_citymodel;

  const GridField2D &_dtm;

  const double top_height;

  const bool fix_buildings;

  // Compute building centroids
  void compute_building_centroids()
  {
    building_centroids.resize(_citymodel.Buildings.size());

    for (size_t i = 0; i < _citymodel.Buildings.size(); i++)
    {
      Vector2D p(0, 0);
      Polygon fp = _citymodel.Buildings[i].Footprint;

      for (auto vertex : fp.Vertices)
      {
        p += Vector2D(vertex);
      }
      p = p / static_cast<double>(fp.Vertices.size());
      building_centroids[i] = p;
    }
  }

  // Compute halos elevation based on min elevation in containing cells
  void compute_halo_elevations()
  {
    std::array<uint, 4> I = {0};

    for (size_t c = 0; c < _mesh.Cells.size(); c++)
    {
      I[0] = _mesh.Cells[c].v0;
      I[1] = _mesh.Cells[c].v1;
      I[2] = _mesh.Cells[c].v2;
      I[3] = _mesh.Cells[c].v3;

      double z_min = std::numeric_limits<double>::max();

      for (size_t i = 0; i < 4; i++)
      {
        const Vector2D p(_mesh.Vertices[I[i]].x, _mesh.Vertices[I[i]].y);
        const double z = _dtm(p);

        z_min = std::min(z_min, z);
      }

      halo_elevations[I[0]] = std::min(halo_elevations[I[0]], z_min);
      halo_elevations[I[1]] = std::min(halo_elevations[I[1]], z_min);
      halo_elevations[I[2]] = std::min(halo_elevations[I[2]], z_min);
      halo_elevations[I[3]] = std::min(halo_elevations[I[3]], z_min);
    }
  }
};

} // namespace DTCC_BUILDER

#endif