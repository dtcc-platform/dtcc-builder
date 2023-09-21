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

#include "StiffnessMatrix.h"
#include "model/GridField.h"

namespace DTCC_BUILDER
{

class BoundaryConditions
{
  typedef unsigned int uint;

public:
  // Vertex Boundary markers:
  // -4 : Neumann vertices
  // -3 : Top Boundary vertices
  // -2 : Ground Boundary vertices
  // -1 : Building Halos Boundary vertices
  std::vector<int> vertex_markers;

  // Boundary values
  std::vector<double> values;

  // Building Polygon Centroids
  std::vector<Vector2D> building_centroids;

  // Elevation for halo vertices based on min Cell elevation
  std::vector<double> halo_elevations;

  // Constructor
  BoundaryConditions(const VolumeMesh &volume_mesh,
                     const City &city,
                     const GridField &dtm,
                     const double top_height,
                     const bool fix_buildings)
      : _volume_mesh(volume_mesh), _city(city), _dtm(dtm),
        top_height(top_height), vertex_markers(volume_mesh.vertices.size(), -4),
        values(volume_mesh.vertices.size(), 0.0), fix_buildings(fix_buildings),
        halo_elevations(volume_mesh.vertices.size(),
                        std::numeric_limits<double>::max())
  {
    // Compute vertex markers
    compute_vertex_markers();

    // Compute boundary values
    compute_boundary_values();
  }

  // Destructor
  ~BoundaryConditions() {}

  // Compute Vertex Boundary markers based on Cell Boundary markers
  void compute_vertex_markers()
  {
    info("Computing vertex markers from cell markers");

    std::array<uint, 4> I = {0};

    // size_t k0 = 0;
    // size_t k1 = 0;
    // size_t k2 = 0;
    // size_t k3 = 0;

    for (size_t c = 0; c < _volume_mesh.cells.size(); c++)
    {
      // Initializing Global Index for each cell
      I[0] = _volume_mesh.cells[c].v0;
      I[1] = _volume_mesh.cells[c].v1;
      I[2] = _volume_mesh.cells[c].v2;
      I[3] = _volume_mesh.cells[c].v3;

      const double z_mean =
          (_volume_mesh.vertices[I[0]].z + _volume_mesh.vertices[I[1]].z +
           _volume_mesh.vertices[I[2]].z + _volume_mesh.vertices[I[3]].z) /
          4;

      const int cell_marker = _volume_mesh.markers[c];
      if (cell_marker >= 0 && fix_buildings) // Building
      {
        for (size_t i = 0; i < 4; i++)
        {
          if (_volume_mesh.vertices[I[i]].z > z_mean)
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
          if (_volume_mesh.vertices[I[i]].z > z_mean)
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
          if (_volume_mesh.vertices[I[i]].z > z_mean)
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
          if (_volume_mesh.vertices[I[i]].z < z_mean)
          {
            continue;
          }
          vertex_markers[I[i]] = std::max(vertex_markers[I[i]], -3);
        }
      }
    }

    // for (size_t v = 0; v < _volume_mesh.vertices.size(); v++)
    // {
    //   if (vertex_markers[v] >= 0)
    //     k0++;
    //   else if (vertex_markers[v] == -1)
    //     k1++;
    //   else if (vertex_markers[v] == -2)
    //     k2++;
    //   else if (vertex_markers[v] == -3)
    //     k3++;
    // }

    return;
  }

  // Compute boundary values
  void compute_boundary_values()
  {
    info("Computing boundary values");

    // TODO: Check if Search tree has already been built
    //_city.build_search_tree(true,0.0);

    // Compute building centroids
    compute_building_centroids();

    // Compute halo elevations
    compute_halo_elevations();

    for (size_t i = 0; i < _volume_mesh.vertices.size(); i++)
    {
      const int vertex_marker = vertex_markers[i];
      if (vertex_marker >= 0) //  && fix_buildings Building
      {
        values[i] = _city.buildings[vertex_marker].max_height() -
                    _volume_mesh.vertices[i].z;
      }
      else if (vertex_marker == -1) // Building Halo
      {
        values[i] = halo_elevations[i] - _volume_mesh.vertices[i].z;
      }
      else if (vertex_marker == -2) // Ground
      {
        const Vector2D p(_volume_mesh.vertices[i].x,
                         _volume_mesh.vertices[i].y);
        values[i] = _dtm(p) - _volume_mesh.vertices[i].z;
      }
      else if (vertex_marker == -3) // Top
      {
        values[i] = top_height - _volume_mesh.vertices[i].z;
      }
      else
      {
        values[i] = 0;
      }
    }
  }

  // Apply boundary conditions on stiffness matrix
  void apply(StiffnessMatrix &A)
  {
    info("Applying boundary conditions to stiffness matrix");

    std::array<uint, 4> I = {0};

    for (size_t c = 0; c < A.shape[0]; c++)
    {
      // Global Index for each cell
      I[0] = _volume_mesh.cells[c].v0;
      I[1] = _volume_mesh.cells[c].v1;
      I[2] = _volume_mesh.cells[c].v2;
      I[3] = _volume_mesh.cells[c].v3;

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

  // Apply boundary conditions on load vector
  void apply(std::vector<double> &b)
  {
    info("Applying boundary conditions to load vector");
    b = values;
  }

private:
  const VolumeMesh &_volume_mesh;

  const City &_city;

  const GridField &_dtm;

  const double top_height;

  const bool fix_buildings;

  // Compute building centroids
  void compute_building_centroids()
  {
    building_centroids.resize(_city.buildings.size());

    for (size_t i = 0; i < _city.buildings.size(); i++)
    {
      Vector2D p(0, 0);
      Polygon fp = _city.buildings[i].footprint;

      for (auto vertex : fp.vertices)
      {
        p += Vector2D(vertex);
      }
      p = p / static_cast<double>(fp.vertices.size());
      building_centroids[i] = p;
    }
  }

  // Compute halos elevation based on min elevation in containing cells
  void compute_halo_elevations()
  {
    std::array<uint, 4> I = {0};

    for (size_t c = 0; c < _volume_mesh.cells.size(); c++)
    {
      I[0] = _volume_mesh.cells[c].v0;
      I[1] = _volume_mesh.cells[c].v1;
      I[2] = _volume_mesh.cells[c].v2;
      I[3] = _volume_mesh.cells[c].v3;

      double z_min = std::numeric_limits<double>::max();

      for (size_t i = 0; i < 4; i++)
      {
        const Vector2D p(_volume_mesh.vertices[I[i]].x,
                         _volume_mesh.vertices[I[i]].y);
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
