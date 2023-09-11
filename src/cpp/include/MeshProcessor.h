// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_PROCESSOR_H
#define DTCC_MESH_PROCESSOR_H

#include <map>
#include <unordered_map>
#include <utility>

#include "Hashing.h"
#include "Logging.h"
#include "model/Mesh.h"
#include "model/VolumeMesh.h"

namespace DTCC_BUILDER
{

class MeshProcessor
{
public:
  /// Compute boundary mesh from volume mesh
  static Mesh compute_boundary_mesh(const VolumeMesh &volume_mesh)
  {
    info("Extracting boundary of 3D mesh...");
    Timer timer("ExtractBoundary3D");

    // Create empty mesh
    Mesh mesh;

    // Map from face --> (num cell neighbors, cell index)
    std::map<Simplex2D, std::pair<size_t, size_t>, CompareSimplex2D> face_map;

    // Iterate over cells and count cell neighbors of faces
    for (size_t i = 0; i < volume_mesh.cells.size(); i++)
    {
      const Simplex3D &cell = volume_mesh.cells[i];
      count_face(face_map, cell.v0, cell.v1, cell.v2, i);
      count_face(face_map, cell.v0, cell.v1, cell.v3, i);
      count_face(face_map, cell.v0, cell.v2, cell.v3, i);
      count_face(face_map, cell.v1, cell.v2, cell.v3, i);
    }

    // Map from old vertex index --> new vertex index
    std::unordered_map<size_t, size_t> vertex_map;

    // Extract faces with exactly one cell neighbor
    for (const auto &it : face_map)
    {
      // Skip if not on boundary
      if (it.second.first != 1)
        continue;

      // Get face and neighboring cell
      const Simplex2D &face = it.first;
      const Simplex3D &cell = volume_mesh.cells[it.second.second];

      // Count vertices (assign new vertex indices)
      const size_t v0 = count_vertex(vertex_map, face.v0);
      const size_t v1 = count_vertex(vertex_map, face.v1);
      const size_t v2 = count_vertex(vertex_map, face.v2);

      // Compute face normal and orientation
      Vector3D n = Geometry::face_normal_3d(face, volume_mesh);
      const Vector3D c = Geometry::cell_center_3d(cell, volume_mesh);
      const Vector3D w = Vector3D(volume_mesh.vertices[face.v0]) - Vector3D(c);
      const int orientation = (Geometry::dot_3d(n, w) > 0.0 ? 1 : -1);

      // Add face and normal
      if (orientation == 1)
      {
        mesh.faces.push_back(Simplex2D{v0, v1, v2});
        mesh.normals.push_back(n);
      }
      else
      {
        mesh.faces.push_back(Simplex2D{v0, v2, v1});
        mesh.normals.push_back(-n);
      }
    }

    // Add vertices
    mesh.vertices.resize(vertex_map.size());
    for (const auto &it : vertex_map)
    {
      const size_t old_index = it.first;
      const size_t new_index = it.second;
      mesh.vertices[new_index] = volume_mesh.vertices[old_index];
    }

    return mesh;
  }

  /// Compute open mesh from boundary, excluding top and sides
  static Mesh compute_open_mesh(const Mesh &boundary)
  {
    info("Comoputing open mesh from boundary...");
    Timer timer("compute_open_mesh");

    // Create empty mesh
    Mesh mesh;

    // Compute bounding box
    BoundingBox3D bbox(boundary.vertices);

    // Map from old vertex index --> new vertex index
    std::unordered_map<size_t, size_t> vertex_map;

    // Extract faces not on top and sides
    for (size_t i = 0; i < boundary.faces.size(); i++)
    {
      // Get face and midpoint
      const Simplex2D &face = boundary.faces[i];
      const Vector3D center = boundary.mid_point(i);

      // Skip if touching bounding box and not pointing upward
      if (std::abs(center.x - bbox.P.x) < Constants::epsilon ||
          std::abs(center.x - bbox.Q.x) < Constants::epsilon ||
          std::abs(center.y - bbox.P.y) < Constants::epsilon ||
          std::abs(center.y - bbox.Q.y) < Constants::epsilon ||
          std::abs(center.z - bbox.Q.z) < Constants::epsilon)
      {
        continue;
      }

      // Count vertices (assign new vertex indices)
      const size_t v0 = count_vertex(vertex_map, face.v0);
      const size_t v1 = count_vertex(vertex_map, face.v1);
      const size_t v2 = count_vertex(vertex_map, face.v2);

      // Add face and normal
      mesh.faces.push_back(Simplex2D{v0, v1, v2});
      mesh.normals.push_back(boundary.normals[i]);
    }

    // Add vertices
    mesh.vertices.resize(vertex_map.size());
    for (const auto &it : vertex_map)
    {
      const size_t old_index = it.first;
      const size_t new_index = it.second;
      mesh.vertices[new_index] = boundary.vertices[old_index];
    }

    return mesh;
  }

  /// Merge meshes into a single mesh
  static Mesh merge_meshes(const std::vector<Mesh> &meshes)
  {
    info("Merging " + str(meshes.size()) + " meshes into a single mesh...");
    Timer timer("merge_meshes");

    // Create empty mesh
    Mesh mesh;

    // Count the number of vertices and cells
    size_t num_vertices = 0;
    size_t num_cells = 0;
    for (size_t i = 0; i < meshes.size(); i++)
    {
      // info("Mesh " + str(i) + " has " + str(meshes[i].vertices.size()) +
      //      " vertices and " + str(meshes[i].faces.size()) + " cells.");
      num_vertices += meshes[i].vertices.size();
      num_cells += meshes[i].faces.size();
    }

    // Allocate arrays
    mesh.vertices.resize(num_vertices);
    mesh.faces.resize(num_cells);

    // Merge data
    size_t k = 0;
    size_t l = 0;
    for (size_t i = 0; i < meshes.size(); i++)
    {
      for (size_t j = 0; j < meshes[i].faces.size(); j++)
      {
        Simplex2D c = meshes[i].faces[j];
        c.v0 += k;
        c.v1 += k;
        c.v2 += k;
        mesh.faces[l++] = c;
      }
      for (size_t j = 0; j < meshes[i].vertices.size(); j++)
        mesh.vertices[k++] = meshes[i].vertices[j];
    }

    return mesh;
  }

private:
  // Count face (number of cell neighbors)
  static void
  count_face(std::map<Simplex2D, std::pair<size_t, size_t>, CompareSimplex2D>
                 &face_map,
             size_t v0,
             size_t v1,
             size_t v2,
             size_t cell_index)
  {
    // Create ordered simplex
    const Simplex2D simplex(v0, v1, v2, true);

    // Check if already added
    auto it = face_map.find(simplex);

    // Add to map or add to counter
    if (it == face_map.end())
      face_map[simplex] = std::make_pair(1, cell_index);
    else if (it->second.first == 1)
      it->second.first++;
    else
    {
      error("Found face with more than 2 cell neighbors.");
    }
  }

  // Count vertex (assign new vertex indices)
  static size_t count_vertex(std::unordered_map<size_t, size_t> &vertex_map,
                             size_t index)
  {
    // Check if already added
    auto it = vertex_map.find(index);

    // Add to map or return existing index
    if (it == vertex_map.end())
    {
      const size_t new_index = vertex_map.size();
      vertex_map[index] = new_index;
      return new_index;
    }
    else
      return it->second;
  }
};

} // namespace DTCC_BUILDER

#endif
