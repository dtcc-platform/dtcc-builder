// Simple vertex-based mesh smoothing (local average)
// Copyright (C) 2020 Anders Logg

#ifndef DTCC_VERTEX_SMOOTHER_H
#define DTCC_VERTEX_SMOOTHER_H

#include <unordered_set>

#include "Surface.h"

namespace DTCC
{

  class VertexSmoother
  {
  public:

    // Smooth 3D surface
    static void SmoothSurface(Surface3D& surface, size_t numSmoothings)
    {
      std::cout << "VertexSmoother: Smoothing surface..." << std::endl;

      // Build vertex connectivity
      std::cout << "VertexSmoother: Building vertex connectivity" << std::endl;
      const size_t numVertices = surface.Vertices.size();
      std::vector<std::unordered_set<size_t>> vertexNeighbors(numVertices);
      for (const auto& T: surface.Cells)
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
        std::cout << "VertexSmoother: Smoothing iteration " << n << std::endl;
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

  };

}

#endif
