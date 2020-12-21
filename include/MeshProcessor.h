// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_PROCESSOR_H
#define DTCC_MESH_PROCESSOR_H

#include <unordered_map>
#include <utility>

#include "Hashing.h"
#include "Logging.h"
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
    surface3D.Faces.clear();
    surface3D.Normals.clear();

    // Map from face --> (num cell neighbors, cell index)
    std::map<Simplex2D, std::pair<size_t, size_t>, CompareSimplex2D> faceMap;

    // Iterate over cells and count cell neighbors of faces
    for (size_t i = 0; i < mesh3D.Cells.size(); i++)
    {
      const Simplex3D &cell = mesh3D.Cells[i];
      CountFace(faceMap, cell.v0, cell.v1, cell.v2, i);
      CountFace(faceMap, cell.v0, cell.v1, cell.v3, i);
      CountFace(faceMap, cell.v0, cell.v2, cell.v3, i);
      CountFace(faceMap, cell.v1, cell.v2, cell.v3, i);
    }

    // Map from old vertex index --> new vertex index
    std::unordered_map<size_t, size_t> vertexMap;

    // Extract faces with exactly one cell neighbor
    for (const auto &it : faceMap)
    {
      // Skip if not on boundary
      if (it.second.first != 1)
        continue;

      // Get face and neighboring cell
      const Simplex2D &face = it.first;
      const Simplex3D &cell = mesh3D.Cells[it.second.second];

      // Count vertices (assign new vertex indices)
      const size_t v0 = CountVertex(vertexMap, face.v0);
      const size_t v1 = CountVertex(vertexMap, face.v1);
      const size_t v2 = CountVertex(vertexMap, face.v2);

      // Compute face normal and orientation
      Vector3D n = Geometry::FaceNormal3D(face, mesh3D);
      const Point3D c = Geometry::CellCenter3D(cell, mesh3D);
      const Vector3D w = Vector3D(mesh3D.Vertices[face.v0]) - Vector3D(c);
      const int orientation = (Geometry::Dot3D(n, w) > 0.0 ? 1 : -1);

      // Add face and normal
      if (orientation == 1)
      {
        surface3D.Faces.push_back(Simplex2D{v0, v1, v2});
        surface3D.Normals.push_back(n);
      }
      else
      {
        surface3D.Faces.push_back(Simplex2D{v0, v2, v1});
        surface3D.Normals.push_back(-n);
      }
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
  static void CountFace(
      std::map<Simplex2D, std::pair<size_t, size_t>, CompareSimplex2D> &faceMap,
      size_t v0,
      size_t v1,
      size_t v2,
      size_t cellIndex)
  {
    // Create ordered simplex
    const Simplex2D simplex(v0, v1, v2, true);

    // Check if already added
    auto it = faceMap.find(simplex);

    // Add to map or add to counter
    if (it == faceMap.end())
      faceMap[simplex] = std::make_pair(1, cellIndex);
    else if (it->second.first == 1)
      it->second.first++;
    else
    {
      Error("Found face with more than 2 cell neighbors.");
    }
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
