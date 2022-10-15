// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_VERTEX_SMOOTHER_H
#define DTCC_VERTEX_SMOOTHER_H

#include <unordered_set>

#include "Logging.h"
#include "Surface.h"
#include "Timer.h"

namespace DTCCBUILDER
{

  class VertexSmoother
  {
  public:
    // Smooth 2D mesh
    static void SmoothMesh(Mesh2D &mesh, size_t numSmoothings)
    {
      info("VertexSmoother: Smoothing mesh...");
      Timer timer("SmoothMesh");

      // Build vertex connectivity
      info("VertexSmoother: Building vertex connectivity");
      const size_t numVertices = mesh.Vertices.size();
      std::vector<std::unordered_set<size_t>> vertexNeighbors(numVertices);
      for (const auto &T : mesh.Cells)
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
        info("VertexSmoother: Smoothing iteration " + str(n));
        for (size_t i = 0; i < numVertices; i++)
        {
          Vector2D p{};
          for (const auto &j : vertexNeighbors[i])
            p += Vector2D(mesh.Vertices[j]);
          p /= static_cast<float>(vertexNeighbors[i].size());
          mesh.Vertices[i] = p;
        }
      }
    }

    // Smooth 3D mesh
    static void SmoothMesh(Mesh3D &mesh, size_t numSmoothings)
    {
      info("VertexSmoother: Smoothing mesh...");
      Timer timer("SmoothMesh");

      // Build vertex connectivity
      info("VertexSmoother: Building vertex connectivity");
      const size_t numVertices = mesh.Vertices.size();
      std::vector<std::unordered_set<size_t>> vertexNeighbors(numVertices);
      for (const auto &T : mesh.Cells)
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
        info("VertexSmoother: Smoothing iteration " + str(n));
        for (size_t i = 0; i < numVertices; i++)
        {
          Vector3D p;
          for (const auto &j : vertexNeighbors[i])
            p += Vector3D(mesh.Vertices[j]);
          p /= static_cast<float>(vertexNeighbors[i].size());
          mesh.Vertices[i] = p;
        }
      }
    }

    // Smooth 3D surface
    static void SmoothSurface(Surface3D& surface, size_t numSmoothings)
    {
      info("VertexSmoother: Smoothing surface...");
      Timer timer("SmoothSurface");

      // Build vertex connectivity
      info("VertexSmoother: Building vertex connectivity");
      const size_t numVertices = surface.Vertices.size();
      std::vector<std::unordered_set<size_t>> vertexNeighbors(numVertices);
      for (const auto &T : surface.Faces)
      {
        vertexNeighbors[T.v0].insert(T.v1);
        vertexNeighbors[T.v0].insert(T.v2);
        vertexNeighbors[T.v1].insert(T.v0);
        vertexNeighbors[T.v1].insert(T.v2);
        vertexNeighbors[T.v2].insert(T.v0);
        vertexNeighbors[T.v2].insert(T.v1);
      }

      // Smooth by setting each vertex coordinate to average of neighbors
      for (size_t n = 0; n < numSmoothings; n++)
      {
        info("VertexSmoother: Smoothing iteration " + str(n));
        for (size_t i = 0; i < numVertices; i++)
        {
          double z = 0.0;
          for (const auto& j: vertexNeighbors[i])
            z += surface.Vertices[j].z;
          z /= static_cast<float>(vertexNeighbors[i].size());
          surface.Vertices[i].z = z;
        }
      }
    }

    // Smooth 2D grid field
    static void SmoothField(GridField2D &field, size_t numSmoothings)
    {
      info("VertexSmoother: Smoothing grid field...");
      Timer timer("SmoothField");

      // Neighbor indices
      std::vector<size_t> indices{};
      indices.reserve(4);

      // Smooth by setting each value to average of neighbors
      for (size_t n = 0; n < numSmoothings; n++)
      {
        info("VertexSmoother: Smoothing iteration " + str(n));
        for (size_t i = 0; i < field.Values.size(); i++)
        {
          // Get neighbors
          indices.clear();
          field.Grid.Index2Boundary(i, indices);

          // Compute average
          double value = 0.0;
          for (const size_t &j : indices)
            value += field.Values[j];
          value /= static_cast<double>(indices.size());
          field.Values[i] = value;
        }
      }
    }
  };

}

#endif
