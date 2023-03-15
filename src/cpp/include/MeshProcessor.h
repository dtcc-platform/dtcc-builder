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

namespace DTCC_BUILDER
{

class MeshProcessor
{
public:
  /// Extract the boundary of a 3D mesh as a 3D surface.
  ///
  /// @param surface The boundary to be extracted as a 3D surface
  /// @param mesh A 3D mesh
  static void ExtractBoundary3D(Surface3D &surface, const Mesh3D &mesh)
  {
    info("MeshProcessor: Extracting boundary of 3D mesh...");
    Timer timer("ExtractBoundary3D");

    // Clear surface
    surface.Vertices.clear();
    surface.Faces.clear();
    surface.Normals.clear();

    // Map from face --> (num cell neighbors, cell index)
    std::map<Simplex2D, std::pair<size_t, size_t>, CompareSimplex2D> faceMap;

    // Iterate over cells and count cell neighbors of faces
    for (size_t i = 0; i < mesh.Cells.size(); i++)
    {
      const Simplex3D &cell = mesh.Cells[i];
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
      const Simplex3D &cell = mesh.Cells[it.second.second];

      // Count vertices (assign new vertex indices)
      const size_t v0 = CountVertex(vertexMap, face.v0);
      const size_t v1 = CountVertex(vertexMap, face.v1);
      const size_t v2 = CountVertex(vertexMap, face.v2);

      // Compute face normal and orientation
      Vector3D n = Geometry::FaceNormal3D(face, mesh);
      const Point3D c = Geometry::CellCenter3D(cell, mesh);
      const Vector3D w = Vector3D(mesh.Vertices[face.v0]) - Vector3D(c);
      const int orientation = (Geometry::Dot3D(n, w) > 0.0 ? 1 : -1);

      // Add face and normal
      if (orientation == 1)
      {
        surface.Faces.push_back(Simplex2D{v0, v1, v2});
        surface.Normals.push_back(n);
      }
      else
      {
        surface.Faces.push_back(Simplex2D{v0, v2, v1});
        surface.Normals.push_back(-n);
      }
    }

    // Add vertices
    surface.Vertices.resize(vertexMap.size());
    for (const auto &it : vertexMap)
    {
      const size_t oldIndex = it.first;
      const size_t newIndex = it.second;
      surface.Vertices[newIndex] = mesh.Vertices[oldIndex];
    }
  }

  /// Extract an open surface from a boundary, excluding top and sides.
  ///
  /// @param surface The open surface to be extracted
  /// @param boundary Closed boundary surface
  static void ExtractOpenSurface3D(Surface3D &surface,
                                   const Surface3D &boundary)
  {
    info("MeshProcessor: Extracting open surface from boundary...");
    Timer timer("ExtractOpenSurface3D");

    // Clear surface
    surface.Vertices.clear();
    surface.Faces.clear();
    surface.Normals.clear();

    // Compute bounding box
    BoundingBox3D bbox(boundary.Vertices);

    // Map from old vertex index --> new vertex index
    std::unordered_map<size_t, size_t> vertexMap;

    // Extract faces not on top and sides
    for (size_t i = 0; i < boundary.Faces.size(); i++)
    {
      // Get face and midpoint
      const Simplex2D &face = boundary.Faces[i];
      const Point3D center = boundary.MidPoint(i);

      // Skip if touching bounding box and not pointing upward
      if (std::abs(center.x - bbox.P.x) < Constants::Epsilon ||
          std::abs(center.x - bbox.Q.x) < Constants::Epsilon ||
          std::abs(center.y - bbox.P.y) < Constants::Epsilon ||
          std::abs(center.y - bbox.Q.y) < Constants::Epsilon ||
          std::abs(center.z - bbox.Q.z) < Constants::Epsilon)
      {
        continue;
      }

      // Count vertices (assign new vertex indices)
      const size_t v0 = CountVertex(vertexMap, face.v0);
      const size_t v1 = CountVertex(vertexMap, face.v1);
      const size_t v2 = CountVertex(vertexMap, face.v2);

      // Add face and normal
      surface.Faces.push_back(Simplex2D{v0, v1, v2});
      surface.Normals.push_back(boundary.Normals[i]);
    }

    // Add vertices
    surface.Vertices.resize(vertexMap.size());
    for (const auto &it : vertexMap)
    {
      const size_t oldIndex = it.first;
      const size_t newIndex = it.second;
      surface.Vertices[newIndex] = boundary.Vertices[oldIndex];
    }
  }

  /// Merge surfaces into a single surface.
  ///
  /// @param surface The single surface
  /// @param mesh3D A
  static void MergeSurfaces3D(Surface3D &surface,
                              const std::vector<Surface3D> &surfaces)
  {
    info("MeshProcessor: Merging 3D surfaces into a single surface...");
    Timer timer("MergeSurfaces3D");

    // Count the number of vertices and cells
    size_t numVertices = 0;
    size_t numCells = 0;
    for (size_t i = 1; i < surfaces.size(); i++)
    {
      numVertices += surfaces[i].Vertices.size();
      numCells += surfaces[i].Faces.size();
    }

    // Allocate arrays
    surface.Vertices.resize(numVertices);
    surface.Faces.resize(numCells);

    // Merge data
    size_t k = 0;
    size_t l = 0;
    for (size_t i = 1; i < surfaces.size(); i++)
    {
      for (size_t j = 0; j < surfaces[i].Faces.size(); j++)
      {
        Simplex2D c = surfaces[i].Faces[j];
        c.v0 += k;
        c.v1 += k;
        c.v2 += k;
        surface.Faces[l++] = c;
      }
      for (size_t j = 0; j < surfaces[i].Vertices.size(); j++)
        surface.Vertices[k++] = surfaces[i].Vertices[j];
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
      error("Found face with more than 2 cell neighbors.");
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

} // namespace DTCC_BUILDER

#endif
