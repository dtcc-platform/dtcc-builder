// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_VERTEX_SMOOTHER_H
#define DTCC_VERTEX_SMOOTHER_H

#include <unordered_set>

#include "Logging.h"
#include "Timer.h"
#include "model/Mesh.h"

namespace DTCC_BUILDER
{

class VertexSmoother
{
public:
  // Smooth mesh
  static void SmoothMesh(Mesh &mesh, size_t numSmoothings)
  {
    info("Smoothing mesh...");
    Timer timer("SmoothMesh");

    // Build vertex connectivity
    info("Building vertex connectivity");
    const size_t numVertices = mesh.Vertices.size();
    std::vector<std::unordered_set<size_t>> vertexNeighbors(numVertices);
    for (const auto &T : mesh.Faces)
    {
      vertexNeighbors[T.v0].insert(T.v1);
      vertexNeighbors[T.v0].insert(T.v2);
      vertexNeighbors[T.v1].insert(T.v2);
      vertexNeighbors[T.v1].insert(T.v0);
      vertexNeighbors[T.v2].insert(T.v0);
      vertexNeighbors[T.v2].insert(T.v1);
    }

    // Smooth by setting each vertex coordinate to average of neighbors
    for (size_t n = 0; n < numSmoothings; n++)
    {
      info("Smoothing iteration " + str(n));
      for (size_t i = 0; i < numVertices; i++)
      {
        Vector3D p{};
        for (const auto &j : vertexNeighbors[i])
          p += Vector3D(mesh.Vertices[j]);
        p /= static_cast<float>(vertexNeighbors[i].size());
        mesh.Vertices[i] = p;
      }
    }
  }

  // Smooth 3D mesh
  static void SmoothMesh(VolumeMesh &volume_mesh, size_t numSmoothings)
  {
    info("Smoothing volume mesh...");
    Timer timer("SmoothMesh");

    // Build vertex connectivity
    info("Building vertex connectivity");
    const size_t numVertices = volume_mesh.Vertices.size();
    std::vector<std::unordered_set<size_t>> vertexNeighbors(numVertices);
    for (const auto &T : volume_mesh.Cells)
    {
      vertexNeighbors[T.v0].insert(T.v1);
      vertexNeighbors[T.v0].insert(T.v2);
      vertexNeighbors[T.v0].insert(T.v3);
      vertexNeighbors[T.v1].insert(T.v2);
      vertexNeighbors[T.v1].insert(T.v3);
      vertexNeighbors[T.v1].insert(T.v0);
      vertexNeighbors[T.v2].insert(T.v3);
      vertexNeighbors[T.v2].insert(T.v0);
      vertexNeighbors[T.v2].insert(T.v1);
      vertexNeighbors[T.v3].insert(T.v0);
      vertexNeighbors[T.v3].insert(T.v1);
      vertexNeighbors[T.v3].insert(T.v2);
    }

    // Smooth by setting each vertex coordinate to average of neighbors
    for (size_t n = 0; n < numSmoothings; n++)
    {
      info("Smoothing iteration " + str(n));
      for (size_t i = 0; i < numVertices; i++)
      {
        Vector3D p;
        for (const auto &j : vertexNeighbors[i])
          p += Vector3D(volume_mesh.Vertices[j]);
        p /= static_cast<float>(vertexNeighbors[i].size());
        volume_mesh.Vertices[i] = p;
      }
    }
  }

  // Smooth grid field
  static GridField smooth_field(const GridField &field, size_t numSmoothings)
  {
    info("Smoothing grid field...");
    Timer timer("smooth_field");

    // Create copy of field
    GridField _field{field};

    // Neighbor indices
    std::vector<size_t> indices{};
    indices.reserve(4);

    // Smooth by setting each value to average of neighbors
    for (size_t n = 0; n < numSmoothings; n++)
    {
      info("Smoothing iteration " + str(n));
      for (size_t i = 0; i < _field.Values.size(); i++)
      {
        // Get neighbors
        indices.clear();
        _field.grid.Index2Boundary(i, indices);

        // Compute average
        double value = 0.0;
        for (const size_t &j : indices)
          value += _field.Values[j];
        value /= static_cast<double>(indices.size());
        _field.Values[i] = value;
      }
    }

    return _field;
  }
};

} // namespace DTCC_BUILDER

#endif
