// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_PROCESSOR_H
#define DTCC_MESH_PROCESSOR_H

#include <unordered_map>
#include <utility>

#include "Hashing.h"
#include "Mesh.h"
#include "Surface.h"

namespace DTCC
{

class MeshProcessor
{
public:
  /// Extract the boundary of a 3D mesh as a 3D surface.
  ///
  /// @param surface3D The boundary to be extracted as a 3D surface
  /// @param mesh3D A 3D mesh
  static void ExtractBoundary3D(Surface3D &surface3D, const Mesh3D &mesh3D)
  {
    Info("MeshProcessor: Extracting boundary of 3D mesh...");
    Timer("ExtractBoundary");

    // Clear surface
    surface3D.Vertices.clear();
    surface3D.Cells.clear();

    // Map from face hash --> face and number of cell neighbors
    std::unordered_map<size_t, std::pair<size_t, Simplex2D>> faceMap;

    // Iterate over cells and count cell neighbors of faces
    for (const auto &cell : mesh3D.Cells)
    {
      CountFace(faceMap, cell.v0, cell.v1, cell.v2);
      CountFace(faceMap, cell.v0, cell.v1, cell.v3);
      CountFace(faceMap, cell.v0, cell.v2, cell.v3);
      CountFace(faceMap, cell.v1, cell.v2, cell.v3);
    }

    // Map from old vertex index --> new vertex index
    std::unordered_map<size_t, size_t> vertexMap;

    // Extract faces with exactly one cell neighbor
    for (const auto &it : faceMap)
    {
      // Skip if not on boundary
      if (it.second.first != 1)
        continue;

      // Get face
      const Simplex2D &simplex = it.second.second;

      // Count vertices (assign new vertex indices)
      const size_t v0 = CountVertex(vertexMap, simplex.v0);
      const size_t v1 = CountVertex(vertexMap, simplex.v1);
      const size_t v2 = CountVertex(vertexMap, simplex.v2);

      // Add face
      surface3D.Cells.push_back(Simplex2D{v0, v1, v2});
    }

    // Add vertices
    surface3D.Vertices.resize(vertexMap.size());
    for (const auto &it : vertexMap)
    {
      const size_t oldIndex = it.first;
      const size_t newIndex = it.second;
      surface3D.Vertices[newIndex] = mesh3D.Vertices[oldIndex];
    }
  }

private:
  // Count face (number of cell neighbors)
  static void
  CountFace(std::unordered_map<size_t, std::pair<size_t, Simplex2D>> &faceMap,
            size_t v0,
            size_t v1,
            size_t v2)
  {
    // Create simplex and hash
    const Simplex2D simplex(v0, v1, v2);
    const size_t hash = Hashing::Hash(simplex);

    // Check if already added
    auto it = faceMap.find(hash);

    // Add to map or add to counter
    if (it == faceMap.end())
      faceMap[hash] = std::make_pair(0, simplex);
    else
      it->second.first++;
  }

  // Count vertex (assign new vertex indices)
  static size_t CountVertex(std::unordered_map<size_t, size_t> &vertexMap,
                            size_t index)
  {
    // Check if already added
    auto it = vertexMap.find(index);

    // Add to map or return existing index
    if (it == vertexMap.end())
    {
      const size_t newIndex = vertexMap.size();
      vertexMap[index] = newIndex;
      return newIndex;
    }
    else
      return it->second;
  }
};

} // namespace DTCC

#endif
