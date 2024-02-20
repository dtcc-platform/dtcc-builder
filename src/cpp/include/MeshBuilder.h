// Copyright (C) 2018 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_BUILDER_H
#define DTCC_MESH_BUILDER_H

#include <cmath>
#include <iostream>
#include <map>
#include <stack>
#include <tuple>
#include <vector>

#include "Eigen/Eigen"
#include "Eigen/Geometry"

#include "BoundingBox.h"
#include "BoundingBoxTree.h"
#include "Geometry.h"
#include "Logging.h"
#include "MeshProcessor.h"
#include "Timer.h"
#include "VertexSmoother.h"
#include "model/GridField.h"
#include "model/Mesh.h"
#include "model/Surface.h"
#include "model/Vector.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern "C"
{
#include <triangle.h>
}

namespace DTCC_BUILDER
{

class MeshBuilder
{
public:
  static Mesh build_terrain_mesh(const std::vector<Polygon> &subdomains,
                                 const GridField &dtm,
                                 double max_mesh_size,
                                 double min_mesh_angle,
                                 size_t smooth_ground = 0)
  {

    // Get bounding box
    const BoundingBox2D &bbox = dtm.grid.bounding_box;
    // build boundary
    Mesh ground_mesh =
        build_ground_mesh(subdomains, bbox.P.x, bbox.P.y, bbox.Q.x, bbox.Q.y,
                          max_mesh_size, min_mesh_angle);
    // Displace ground surface. Fill all points with maximum height. This is
    // used to always choose the smallest height for each point since each point
    // may be visited multiple times.
    const double z_max = dtm.max();
    for (size_t i = 0; i < ground_mesh.vertices.size(); i++)
      ground_mesh.vertices[i].z = z_max;

    // If ground is not float, iterate over the triangles
    for (size_t i = 0; i < ground_mesh.faces.size(); i++)
    {
      // Get cell marker
      const int cell_marker = ground_mesh.markers[i];

      // Get triangle
      const Simplex2D &T = ground_mesh.faces[i];

      // Check cell marker
      if (cell_marker != -2) // not ground
      {
        // Compute minimum height of vertices
        double z_min = std::numeric_limits<double>::max();
        z_min = std::min(z_min, dtm(ground_mesh.vertices[T.v0]));
        z_min = std::min(z_min, dtm(ground_mesh.vertices[T.v1]));
        z_min = std::min(z_min, dtm(ground_mesh.vertices[T.v2]));

        // Set minimum height for all vertices
        set_min(ground_mesh.vertices[T.v0].z, z_min);
        set_min(ground_mesh.vertices[T.v1].z, z_min);
        set_min(ground_mesh.vertices[T.v2].z, z_min);
      }
      else
      {
        // Sample height map at vertex position for all vertices
        set_min(ground_mesh.vertices[T.v0].z, dtm(ground_mesh.vertices[T.v0]));
        set_min(ground_mesh.vertices[T.v1].z, dtm(ground_mesh.vertices[T.v1]));
        set_min(ground_mesh.vertices[T.v2].z, dtm(ground_mesh.vertices[T.v2]));
      }
    }
    info("smooth ground...");
    if (smooth_ground > 0)
      VertexSmoother::smooth_mesh(ground_mesh, smooth_ground, true);

    // for (size_t i = 0; i < ground_mesh.faces.size(); i++)
    // {
    //   auto normal = Geometry::face_normal(ground_mesh.faces[i], ground_mesh);
    //   if (normal.z < 0)
    //     ground_mesh.faces[i].flip();
    // }

    info("ground mesh done");
    return ground_mesh;
  }

  // build ground mesh for city.
  //
  // The mesh is a triangular mesh of the rectangular region
  // defined by (xmin, xmax) x (ymin, ymax). The edges of the mesh respect
  // the boundaries of the buildings.
  //
  // markers:
  //
  // -2: ground (cells outside buildings and halos)
  // -1: halos (cells close to buildings)
  //  0: building 0 (cells inside building 0)
  //  1: building 1 (cells inside building 1)
  //  etc (non-negative integers mark cells inside buildings)
  static Mesh build_ground_mesh(const std::vector<Polygon> &subdomains,
                                double xmin,
                                double ymin,
                                double xmax,
                                double ymax,
                                double max_mesh_size,
                                double min_mesh_angle)
  {
    info("Building ground mesh for city...");
    Timer timer("build_ground_mesh");

    // print some stats
    const BoundingBox2D bounding_box(Vector2D(xmin, ymin),
                                     Vector2D(xmax, ymax));
    const size_t nx = (bounding_box.Q.x - bounding_box.P.x) / max_mesh_size;
    const size_t ny = (bounding_box.Q.y - bounding_box.P.y) / max_mesh_size;
    const size_t n = nx * ny;
    info("Domain bounding box is " + str(bounding_box));
    info("Maximum mesh size is " + str(max_mesh_size));
    info("Estimated number of triangles is " + str(n));
    info("Number of subdomains (buildings) is " + str(subdomains.size()));

    // Extract subdomains (building footprints)
    std::vector<std::vector<Vector2D>> triangle_sub_domains;
    for (auto const &sd : subdomains)
    {
      triangle_sub_domains.push_back(sd.vertices);
      for (auto const &hole : sd.holes)
        triangle_sub_domains.push_back(hole);
    }
    info("Number of subdomains (buildings + holes) is " +
         str(triangle_sub_domains.size()));
    // build boundary
    std::vector<Vector2D> boundary{};
    boundary.push_back(bounding_box.P);
    boundary.push_back(Vector2D(bounding_box.Q.x, bounding_box.P.y));
    boundary.push_back(bounding_box.Q);
    boundary.push_back(Vector2D(bounding_box.P.x, bounding_box.Q.y));

    // build 2D mesh
    Mesh mesh;
    call_triangle(mesh, boundary, triangle_sub_domains, max_mesh_size,
                  min_mesh_angle, false);

    // Mark subdomains
    compute_domain_markers(mesh, subdomains);

    return mesh;
  }

  // Layer ground mesh to create a volume mesh.
  //
  // The volume mesh is a tetrahedral mesh constructed
  // extruding the 2D mesh in the vertical (z) direction.
  //
  // markers:
  //
  // -2: ground (cells outside buildings and halos)
  // -1: halos (cells close to buildings)
  //  0: building 0 (cells inside building 0)
  //  1: building 1 (cells inside building 1)
  //  etc (non-negative integers mark cells inside buildings)
  //
  // Note that the markers are just propagatated upward from the
  // ground mesh, meaning that the markers will be the same in each
  // column of the ground mesh.
  //
  // It is assumed that the ground mesh is sorted, which is the case if the
  // mesh has been built by calling build_mesh().
  static VolumeMesh layer_ground_mesh(const Mesh &ground_mesh,
                                      double domain_height,
                                      double max_mesh_size)
  {
    Timer timer("build_volume_mesh");

    // Compute number of layers
    const size_t num_layers = int(std::ceil(domain_height / max_mesh_size));
    const double dz = domain_height / double(num_layers);
    const size_t layer_size = ground_mesh.vertices.size();

    info("Building volume mesh with " + str(num_layers) + " layers...");

    // Initialize volume mesh
    VolumeMesh volume_mesh;
    volume_mesh.vertices.resize((num_layers + 1) * ground_mesh.vertices.size());
    volume_mesh.cells.resize(num_layers * 3 * ground_mesh.faces.size());
    volume_mesh.markers.resize(volume_mesh.cells.size());
    volume_mesh.num_layers = num_layers;

    // Add vertices
    {
      size_t k = 0;
      for (size_t layer = 0; layer <= num_layers; layer++)
      {
        // Compute height of layer
        const double z = layer * dz;

        // Iterate over vertices in layer
        for (const auto &p_2d : ground_mesh.vertices)
          volume_mesh.vertices[k++] = Vector3D(p_2d.x, p_2d.y, z);
      }
    }

    // Add cells
    {
      size_t k = 0;
      size_t offset = 0;
      for (size_t layer = 0; layer < num_layers; layer++)
      {
        // Iterate over triangles in layer
        for (const auto &T : ground_mesh.faces)
        {
          // Get sorted vertex indices for bottom layer
          const size_t u0 = T.v0 + offset;
          const size_t u1 = T.v1 + offset;
          const size_t u2 = T.v2 + offset;

          // Get sorted vertices for top layer
          const size_t v0 = u0 + layer_size;
          const size_t v1 = u1 + layer_size;
          const size_t v2 = u2 + layer_size;

          // Create three tetrahedra by connecting the first vertex
          // of each edge in the bottom layer with the second
          // vertex of the corresponding edge in the top layer.
          volume_mesh.cells[k++] = Simplex3D(u0, u1, u2, v2);
          volume_mesh.cells[k++] = Simplex3D(u0, v1, u1, v2);
          volume_mesh.cells[k++] = Simplex3D(u0, v0, v1, v2);
        }

        // Add to offset
        offset += layer_size;
      }
    }

    // Add domain markers
    {
      size_t k = 0;
      for (size_t layer = 0; layer < num_layers; layer++)
      {
        for (const auto &marker : ground_mesh.markers)
        {
          int m = 0;

          // Top layer marked as -3
          if (layer == num_layers - 1)
            m = -3;

          // Halo and ground only marked for bottom layer
          else if (marker == -1 || marker == -2)
          {
            if (layer == 0)
              m = marker;
            else
              m = -4;
          }

          // buildings marked for all layers (except top layer).
          // Later adjusted to -4 above buildings in trim_volume_mesh.
          else
          {
            m = marker;
          }

          // Set markers
          volume_mesh.markers[k++] = m;
          volume_mesh.markers[k++] = m;
          volume_mesh.markers[k++] = m;
        }
      }
    }

    return volume_mesh;
  }

  // Trim volume mesh by removing cells inside buildings.
  //
  // markers:
  //
  // -4: other (not top, ground, halo, or building)
  // -3: top (cells in *top layer*)
  // -2: ground (cells in *bottom layer* outside buildings and halos)
  // -1: halos (cells in *bottom layer* close to buildings)
  //  0: building 0 (cells *first layer* inside building 0)
  //  1: building 1 (cells *first layer* inside building 1)
  //  etc (non-negative integers mark cells inside buildings)
  //
  // Note that the markers are adjusted in the following way:
  //
  // - Only cells in bottom layer marked as ground (-2) or halo (-1)
  // - Only cells in first layer above a building marked as building
  // - cells in top layer marked as top (-3)
  // - All other cells (in between) marked as other (-4)
  static VolumeMesh trim_volume_mesh(const VolumeMesh &volume_mesh,
                                     const Mesh &mesh,
                                     const std::vector<Surface> &buildings)
  {
    info("Trimming volume mesh...");
    Timer timer("trim_volume_mesh");

    // Get sizes
    const size_t num_buildings = buildings.size();
    const size_t num_cells_2d = mesh.faces.size();
    const size_t num_cells_3d = volume_mesh.cells.size();
    const size_t layer_size = 3 * mesh.faces.size();

    // Phase 1: Determine which cells should be trimmed
    // ------------------------------------------------

    // build map from buildings to cells in 2D mesh
    std::vector<std::vector<size_t>> building_cells_2d(num_buildings);
    for (size_t cell_index_2d = 0; cell_index_2d < num_cells_2d;
         cell_index_2d++)
    {
      const int building_index = mesh.markers[cell_index_2d];
      if (building_index >= 0)
        building_cells_2d[building_index].push_back(cell_index_2d);
    }

    // Create markers for cells to be trimmed (keep by default)
    std::vector<bool> trim_cell(num_cells_3d);
    std::fill(trim_cell.begin(), trim_cell.end(), false);

    // Keep track of first layer for each building
    std::vector<size_t> first_layer(num_buildings);
    std::fill(first_layer.begin(), first_layer.end(), 0);

    // Iterate over buildings
    for (size_t building_index = 0; building_index < num_buildings;
         building_index++)
    {
      // Iterate over layers
      for (size_t layer = 0; layer < volume_mesh.num_layers; layer++)
      {
        // build list of 3D cells for building in current layer
        std::vector<size_t> cells_3d;
        for (const auto &cell_index_2d : building_cells_2d[building_index])
        {
          for (size_t j = 0; j < 3; j++)
            cells_3d.push_back(index_3d(layer, layer_size, cell_index_2d, j));
        }

        // Trim layer if any cell midpoint is below building height
        bool trim_layer = false;
        for (const auto &cell_index_3d : cells_3d)
        {
          const double z = volume_mesh.mid_point(cell_index_3d).z;
          const double h = buildings[building_index].max_height();
          if (z < h)
          {
            trim_layer = true;
            break;
          }
        }

        // Check if layer should be trimmed
        if (trim_layer)
        {
          // Mark cells for trimming
          for (const auto &cell_index_3d : cells_3d)
            trim_cell[cell_index_3d] = true;
        }
        else
        {
          // If layer should be kept, no need to check more layers
          first_layer[building_index] = layer;
          break;
        }
      }
    }

    // Phase 2: Adjust markers
    // -----------------------

    // Create copy of markeres
    std::vector<int> markers{volume_mesh.markers};

    // Mark cells between bottom and top layer as -4
    for (size_t layer = 1; layer < volume_mesh.num_layers - 1; layer++)
    {
      for (size_t cell_index_2d = 0; cell_index_2d < num_cells_2d;
           cell_index_2d++)
        for (size_t j = 0; j < 3; j++)
          markers[index_3d(layer, layer_size, cell_index_2d, j)] = -4;
    }

    // Mark cells in top layer as -3
    for (size_t cell_index_2d = 0; cell_index_2d < num_cells_2d;
         cell_index_2d++)
    {
      for (size_t j = 0; j < 3; j++)
        markers[index_3d(volume_mesh.num_layers - 1, layer_size, cell_index_2d,
                         j)] = -3;
    }

    // Mark cells in first layer above each building:
    //
    // 0, 1, 2, ... if building is not covered by bottom layer (normal case)
    // -1           if building is covered by bottom layer (modify to halo)
    for (size_t building_index = 0; building_index < num_buildings;
         building_index++)
    {
      const size_t layer = first_layer[building_index];
      size_t marker = building_index;
      if (layer == 0)
      {
        warning("Building " + str(building_index) +
                " is covered by bottom layer");
        marker = -1;
      }
      for (const auto &cell_index_2d : building_cells_2d[building_index])
      {
        for (size_t j = 0; j < 3; j++)
          markers[index_3d(layer, layer_size, cell_index_2d, j)] = marker;
      }
    }

    // Phase 3: Extract new mesh for all cells that should be kept
    // -----------------------------------------------------------

    // Renumber vertices and cells
    std::unordered_map<size_t, size_t> vertex_map;
    std::unordered_map<size_t, size_t> cell_map;
    size_t k = 0;
    size_t l = 0;
    for (size_t cell_index_3d = 0; cell_index_3d < volume_mesh.cells.size();
         cell_index_3d++)
    {
      if (!trim_cell[cell_index_3d])
      {
        // Get cell
        const Simplex3D &T = volume_mesh.cells[cell_index_3d];

        // Renumbers vertices
        if (vertex_map.find(T.v0) == vertex_map.end())
          vertex_map[T.v0] = k++;
        if (vertex_map.find(T.v1) == vertex_map.end())
          vertex_map[T.v1] = k++;
        if (vertex_map.find(T.v2) == vertex_map.end())
          vertex_map[T.v2] = k++;
        if (vertex_map.find(T.v3) == vertex_map.end())
          vertex_map[T.v3] = k++;

        // Renumber cells
        cell_map[cell_index_3d] = l++;
      }
    }

    // Initialize new mesh data
    const size_t num_vertices = vertex_map.size();
    const size_t num_cells = cell_map.size();
    std::vector<Vector3D> _vertices(num_vertices);
    std::vector<Simplex3D> _cells(num_cells);
    std::vector<int> _markers(num_cells);

    // Set new mesh data
    for (const auto v : vertex_map)
      _vertices[v.second] = volume_mesh.vertices[v.first];
    for (const auto c : cell_map)
    {
      _cells[c.second].v0 = vertex_map[volume_mesh.cells[c.first].v0];
      _cells[c.second].v1 = vertex_map[volume_mesh.cells[c.first].v1];
      _cells[c.second].v2 = vertex_map[volume_mesh.cells[c.first].v2];
      _cells[c.second].v3 = vertex_map[volume_mesh.cells[c.first].v3];
      _markers[c.second] = markers[c.first];
    }

    // Create new mesh and assign data
    VolumeMesh _volume_mesh;
    _volume_mesh.vertices = _vertices;
    _volume_mesh.cells = _cells;
    _volume_mesh.markers = _markers;

    return _volume_mesh;
  }

  static std::vector<Mesh>
  build_city_surface_mesh(const std::vector<Surface> &buildings,
                          const GridField &dtm,
                          double max_mesh_size,
                          double min_mesh_angle,
                          size_t smooth_ground = 0,
                          bool merge_meshes = true)
  {
    auto build_city_surface_t = Timer("build_city_surface_mesh");
    auto terrain_time = Timer("build_city_surface_mesh: step 1 terrain");
    std::vector<Polygon> subdomains;
    for (const auto &b : buildings)
    {
      subdomains.push_back(b.to_polygon());
    }
    Mesh terrain_mesh = build_terrain_mesh(subdomains, dtm, max_mesh_size,
                                           min_mesh_angle, smooth_ground);
    terrain_time.stop();

    std::vector<Mesh> city_mesh;
    std::vector<Mesh> building_meshes;

    std::map<size_t, std::vector<Simplex2D>> building_faces;
    std::vector<size_t> building_indices;
    info("finding markes");
    auto find_markers_t = Timer("build_city_surface_mesh: step 2 find markers");
    for (size_t i = 0; i < terrain_mesh.markers.size(); i++)
    {
      auto marker = terrain_mesh.markers[i];

      if (marker >= 0)
      {
        // info("marker: " + str(marker) + " i: " + str(i));
        building_faces[marker].push_back(terrain_mesh.faces[i]);
        building_indices.push_back(i);
      }
    }
    find_markers_t.stop();

    info("building meshes");
    auto building_meshes_t =
        Timer("build_city_surface_mesh: step 3 building meshes");
    for (auto it = building_faces.begin(); it != building_faces.end(); ++it)
    {
      auto marker = it->first;
      auto faces = it->second;
      auto building = buildings[marker];
      auto roof_height = building.max_height();

      Mesh building_mesh;

      // add roofs
      for (const auto &face : faces)
      {
        auto v0 = terrain_mesh.vertices[face.v0];
        auto v1 = terrain_mesh.vertices[face.v1];
        auto v2 = terrain_mesh.vertices[face.v2];
        auto v3 = Vector3D(v0.x, v0.y, roof_height);
        auto v4 = Vector3D(v1.x, v1.y, roof_height);
        auto v5 = Vector3D(v2.x, v2.y, roof_height);

        auto num_vertices = building_mesh.vertices.size();
        building_mesh.vertices.push_back(v3);
        building_mesh.vertices.push_back(v4);
        building_mesh.vertices.push_back(v5);
        // info("vertices: " + str(v3) + " " + str(v4) + " " + str(v5));
        // info("faces: " + str(face.v0) + " " + str(face.v1) + " " +
        //     str(face.v2));
        building_mesh.faces.push_back(
            Simplex2D(num_vertices, num_vertices + 1, num_vertices + 2));
      }

      // add walls
      auto naked_edges = MeshProcessor::find_naked_edges(faces);
      for (const auto &edge_faces : naked_edges)
      {
        Simplex1D edge = edge_faces.first;
        Simplex2D edge_face = edge_faces.second;
        auto face_center = Geometry::face_center(edge_face, terrain_mesh);

        auto ground_v0 = terrain_mesh.vertices[edge.v0];
        auto ground_v1 = terrain_mesh.vertices[edge.v1];
        auto roof_v0 = Vector3D(ground_v0.x, ground_v0.y, roof_height);
        auto roof_v1 = Vector3D(ground_v1.x, ground_v1.y, roof_height);

        Mesh wall_mesh;
        wall_mesh.vertices = {ground_v0, ground_v1, roof_v0, roof_v1};
        auto wall_normal =
            Geometry::triangle_normal(ground_v0, ground_v1, roof_v1);
        if (Geometry::dot_3d(wall_normal, face_center - ground_v0) > 0)
        {
          wall_mesh.faces = {Simplex2D(0, 1, 3), Simplex2D(0, 3, 2)};
        }
        else
        {
          wall_mesh.faces = {Simplex2D(0, 3, 1), Simplex2D(0, 2, 3)};
        }

        building_mesh = MeshProcessor::merge_meshes({building_mesh, wall_mesh});
      }

      building_meshes.push_back(MeshProcessor::weld_mesh(building_mesh));
    }
    building_meshes_t.stop();

    // remove triangles inside houses from terrain
    auto remove_inside_t =
        Timer("build_city_surface_mesh: step 4 remove inside");
    std::sort(building_indices.begin(), building_indices.end());

    // if index is in list of building indices, move to the end of the list
    auto new_end = std::remove_if(
        terrain_mesh.faces.begin(), terrain_mesh.faces.end(),
        [&](const auto &face)
        {
          return std::binary_search(
              building_indices.begin(), building_indices.end(),
              &face -
                  &terrain_mesh.faces[0]); // calculate the index of the face in
                                           // the terrain_mesh.faces vector.
        });
    // remove all elements that have been moved
    terrain_mesh.faces.erase(new_end, terrain_mesh.faces.end());
    remove_inside_t.stop();

    auto final_merger_t = Timer("build_city_surface_mesh: step 5 final merge");
    city_mesh.push_back(terrain_mesh);
    city_mesh.insert(city_mesh.end(), building_meshes.begin(),
                     building_meshes.end());
    if (merge_meshes)
    {
      auto merged_mesh = MeshProcessor::merge_meshes(city_mesh, true);
      city_mesh = {merged_mesh};
    }
    final_merger_t.stop();
    build_city_surface_t.stop();
    // Timer::report("city surface");
    return city_mesh;
  }

  static Mesh mesh_surface(const Surface &surface,

                           double max_triangle_area_size = -1,
                           double min_mesh_angle = 25)
  // convert 3D Surface to triangle Mesh. If max_triangle_area_size is
  // greater than 0, triangle will be used for triangulations, otherwise
  // it will be naively triangulated by connecting the vertices in order.
  {
    Mesh mesh;
    if (surface.vertices.size() < 3)
      return mesh;
    if (max_triangle_area_size <= 0)
    {
      fast_mesh(mesh, surface);
    }
    else
    {
      call_triangle(mesh, surface, max_triangle_area_size, min_mesh_angle);
    }
    return mesh;
  }

  static Mesh mesh_multisurface(const MultiSurface &multi_surface,
                                double max_triangle_area_size = -1,
                                double min_mesh_angle = 25,
                                bool weld = false)
  {
    std::vector<Mesh> multimesh(multi_surface.surfaces.size());
#pragma omp parallel for
    for (size_t i = 0; i < multi_surface.surfaces.size(); i++)
    {
      multimesh[i] = mesh_surface(multi_surface.surfaces[i],
                                  max_triangle_area_size, min_mesh_angle);
    }
    // for (const auto &surface : multi_surface.surfaces)
    // {
    //   auto surface_mesh =
    //       mesh_surface(surface, max_triangle_area_size, min_mesh_angle);
    //   multimesh.push_back(surface_mesh);
    // }
    auto mesh = MeshProcessor::merge_meshes(multimesh, weld);
    return mesh;
  }

  static std::vector<Mesh>
  mesh_multisurfaces(const std::vector<MultiSurface> &multi_surfaces,
                     double max_triangle_area_size = -1,
                     double min_mesh_angle = 25,
                     bool weld = false)
  {
    size_t n = multi_surfaces.size();
    std::vector<Mesh> meshes(n);
#pragma omp parallel for
    for (size_t i = 0; i < n; i++)
    {
      auto mesh = mesh_multisurface(multi_surfaces[i], max_triangle_area_size,
                                    min_mesh_angle, weld);
      meshes[i] = mesh;
    }

    return meshes;
  }

private:
  // Map from 2D cell index to 3D cell indices
  static size_t
  index_3d(size_t layer, size_t layer_size, size_t cell_index_2d, size_t j)
  {
    return layer * layer_size + 3 * cell_index_2d + j;
  }

  static void call_triangle(Mesh &mesh,
                            const Surface &surface,
                            double max_mesh_size,
                            double min_mesh_angle)
  {
    std::vector<std::vector<Vector2D>> sd;

    const auto z_normal = Eigen::Vector3d(0, 0, 1);
    auto normal = Geometry::surface_noraml(surface);
    auto e_norm = Eigen::Vector3d(normal.x, normal.y, normal.z);
    auto centroid = Geometry::surface_centroid(surface);
    auto e_centroid = Eigen::Vector3d(centroid.x, centroid.y, centroid.z);

    // auto rot_matrix =
    //     Eigen::Matrix3d(Eigen::Quaterniond::FromTwoVectors(e_norm,
    //     z_normal));
    // std::cout << "rot_matrix: " << rot_matrix << std::endl;
    auto transform = Eigen::Transform<double, 3, Eigen::Isometry>();
    auto translation = Eigen::Translation3d(-e_centroid);
    auto rotation = Eigen::Quaterniond::FromTwoVectors(e_norm, z_normal);

    transform = rotation * translation;
    auto transform_inv = transform.inverse();
    // std::cout << "trans_matrix " << trans_matrix << std::endl;

    // std::cout << "transform " << transform << std::endl;
    // auto transform_inv = transform.inverse();
    // std::cout << "transform inv" << transform_inv << std::endl;
    std::vector<Vector2D> projected_surface;
    for (const auto &v : surface.vertices)
    {
      auto e_v = Eigen::Vector3d(v.x, v.y, v.z);
      auto e_v_prime = transform * e_v;
      projected_surface.push_back(Vector2D(e_v_prime.x(), e_v_prime.y()));
    }
    call_triangle(mesh, projected_surface, sd, max_mesh_size, min_mesh_angle,
                  false);
    for (auto &v : mesh.vertices)
    {
      auto e_v = Eigen::Vector3d(v.x, v.y, 0);
      auto e_v_prime = transform_inv * e_v;
      v.x = e_v_prime.x();
      v.y = e_v_prime.y();
      v.z = e_v_prime.z();
    }
  }

  static void fast_mesh(Mesh &mesh, const Surface &surface)
  {
    if (surface.vertices.size() < 3)
      return;
    if (surface.vertices.size() == 3)
    {
      mesh.vertices = surface.vertices;
      mesh.faces.push_back(Simplex2D(0, 1, 2));
      return;
    }

    if (Geometry::is_convex(surface))
    {
      mesh.vertices = surface.vertices;
      for (size_t i = 1; i < surface.vertices.size() - 1; i++)
      {
        mesh.faces.push_back(Simplex2D(0, i, i + 1));
      }
    }
    else
    {
      call_triangle(mesh, surface, -1, -1);
    }
  }

  // Call Triangle to compute 2D mesh
  static void
  call_triangle(Mesh &mesh,
                const std::vector<Vector2D> &boundary,
                const std::vector<std::vector<Vector2D>> &sub_domains,
                double max_mesh_size,
                double min_mesh_angle,
                bool sort_triangles)
  {
    Timer timer("call_triangle");

    // Set area constraint to control mesh size
    const double max_area = 0.5 * max_mesh_size * max_mesh_size;

    // Set input switches for Triangle
    std::string triswitches = "zQp";
    if (min_mesh_angle > 0)
      triswitches += "q" + str(min_mesh_angle, 3);
    if (max_mesh_size > 0)
      triswitches += "a" + str(max_area, 3);
    debug("Triangle switches: " + triswitches);

    // Convert to C-style string
    char *triswitches_c = new char[triswitches.length() + 1];
    std::strcpy(triswitches_c, triswitches.c_str());

    // z = use zero-based numbering
    // p = use polygon input (segments)
    // q = control mesh quality
    //
    // Note that the minimum angle (here 25) should be
    // as large as possible for high quality meshes but
    // it should be less than 28.6 degrees to guarantee
    // that Triangle terminates. Default is 20 degrees.

    // Create input data structure for Triangle
    struct triangulateio in = create_triangle_io();

    // Set number of points
    size_t num_points = boundary.size();
    for (auto const &innerPolygon : sub_domains)
    {
      num_points += innerPolygon.size();
    }
    in.numberofpoints = num_points;

    // Set points
    in.pointlist = new double[2 * num_points];
    {
      size_t k = 0;
      for (auto const &p : boundary)
      {
        in.pointlist[k++] = p.x;
        in.pointlist[k++] = p.y;
      }
      for (auto const &innerPolygon : sub_domains)
      {
        for (auto const &p : innerPolygon)
        {
          in.pointlist[k++] = p.x;
          in.pointlist[k++] = p.y;
        }
      }
    }

    // Set number of segments
    const size_t num_segments = num_points;
    in.numberofsegments = num_segments;

    // Set segments
    in.segmentlist = new int[2 * num_segments];
    {
      size_t k = 0;
      size_t n = 0;
      for (size_t j = 0; j < boundary.size(); j++)
      {
        const size_t j0 = j;
        const size_t j1 = (j + 1) % boundary.size();
        in.segmentlist[k++] = n + j0;
        in.segmentlist[k++] = n + j1;
      }
      n += boundary.size();
      for (size_t i = 0; i < sub_domains.size(); i++)
      {
        for (size_t j = 0; j < sub_domains[i].size(); j++)
        {
          const size_t j0 = j;
          const size_t j1 = (j + 1) % sub_domains[i].size();
          in.segmentlist[k++] = n + j0;
          in.segmentlist[k++] = n + j1;
        }
        n += sub_domains[i].size();
      }
    }

    // Note: This is how set holes but it's not used here since we
    // need the triangles for the interior *above* the buildings.

    /*
    // Set number of holes
    const size_t numHoles = SubDomains.size();
    in.numberofholes = numHoles;

    // Set holes. Note that we assume that we can get an
    // interior point of each hole (inner polygon) by computing
    // its center of mass.
    in.holelist = new double[2 * numHoles];
    {
    size_t k = 0;
    Vector2D c;
    for (auto const & InnerPolygon : SubDomains)
    {
    for (auto const & p : InnerPolygon)
    {
    c += p;
    }
    c /= InnerPolygon.size();
    in.holelist[k++] = c.x;
    in.holelist[k++] = c.y;
    }
    }
    */

    // Prepare output data for Triangl;e
    struct triangulateio out = create_triangle_io();
    struct triangulateio vorout = create_triangle_io();

    // Call Triangle
    triangulate(triswitches_c, &in, &out, &vorout);
    delete[] triswitches_c;

    // Uncomment for debugging
    // print_triangle_io(out);
    // print_triangle_io(vorout);

    // Extract points
    mesh.vertices.reserve(out.numberofpoints);
    for (int i = 0; i < out.numberofpoints; i++)
    {
      Vector3D p(out.pointlist[2 * i], out.pointlist[2 * i + 1], 0.0);
      mesh.vertices.push_back(p);
    }

    // Extract triangles
    mesh.faces.reserve(out.numberoftriangles);
    for (int i = 0; i < out.numberoftriangles; i++)
    {
      // Note the importance of creating a sorted simplex here!
      Simplex2D t(out.trianglelist[3 * i], out.trianglelist[3 * i + 1],
                  out.trianglelist[3 * i + 2], sort_triangles);
      mesh.faces.push_back(t);
    }

    // Free memory
    // trifree(&out); // causes segfault
    delete[] in.pointlist;
    delete[] in.segmentlist;
    delete[] in.holelist;
  }

  // Create and reset Triangle I/O data structure
  static struct triangulateio create_triangle_io()
  {
    struct triangulateio io;

    io.pointlist = nullptr;
    io.pointmarkerlist = nullptr;
    io.pointmarkerlist = nullptr;
    io.numberofpoints = 0;
    io.numberofpointattributes = 0;
    io.trianglelist = nullptr;
    io.triangleattributelist = nullptr;
    io.trianglearealist = nullptr;
    io.neighborlist = nullptr;
    io.numberoftriangles = 0;
    io.numberofcorners = 0;
    io.numberoftriangleattributes = 0;
    io.segmentlist = nullptr;
    io.segmentmarkerlist = nullptr;
    io.numberofsegments = 0;
    io.holelist = nullptr;
    io.numberofholes = 0;
    io.regionlist = nullptr;
    io.numberofregions = 0;
    io.edgelist = nullptr;
    io.edgemarkerlist = nullptr;
    io.normlist = nullptr;
    io.numberofedges = 0;

    return io;
  }

  // print triangle I/O data
  static void print_triangle_io(const struct triangulateio &io)
  {
    info("Triangle I/O data: ");
    info("  pointlist = " +
         str(reinterpret_cast<std::uintptr_t>(io.pointlist)));
    info("  pointmarkerlist = " +
         str(reinterpret_cast<std::uintptr_t>(io.pointmarkerlist)));
    if (io.pointmarkerlist)
    {
      std::stringstream string_builder{};
      string_builder << "   ";
      for (int i = 0; i < io.numberofpoints; i++)
        string_builder << " " << io.pointmarkerlist[i];
      string_builder << std::endl;
      info(string_builder.str());
    }
    info("  numberofpoints = " + str(io.numberofpoints));
    info("  numberofpointattributes = " + str(io.numberofpointattributes));
    info("  trianglelist = " +
         str(reinterpret_cast<std::uintptr_t>(io.trianglelist)));
    info("  triangleattributelist = " +
         str(reinterpret_cast<std::uintptr_t>(io.triangleattributelist)));
    info("  trianglearealist = " +
         str(reinterpret_cast<std::uintptr_t>(io.trianglearealist)));
    info("  neighborlist = " +
         str(reinterpret_cast<std::uintptr_t>(io.neighborlist)));
    info("  numberoftriangles = " + str(io.numberoftriangles));
    info("  numberofcorners = " + str(io.numberofcorners));
    info("  numberoftriangleattributes = " +
         str(io.numberoftriangleattributes));
    info("  segmentlist = " +
         str(reinterpret_cast<std::uintptr_t>(io.segmentlist)));
    info("  segmentmarkerlist = " +
         str(reinterpret_cast<std::uintptr_t>(io.segmentmarkerlist)));
    if (io.segmentmarkerlist)
    {
      std::stringstream string_builder{};
      string_builder << "   ";
      for (int i = 0; i < io.numberofsegments; i++)
        string_builder << " " << io.segmentmarkerlist[i];
      string_builder << std::endl;
      info(string_builder.str());
    }
    info("  numberofsegments = " + str(io.numberofsegments));
    info("  holelist = " + str(reinterpret_cast<std::uintptr_t>(io.holelist)));
    info("  numberofholes = " + str(io.numberofholes));
    info("  regionlist = " +
         str(reinterpret_cast<std::uintptr_t>(io.regionlist)));
    info("  numberofregions = " + str(io.numberofregions));
    info("  edgelist = " + str(reinterpret_cast<std::uintptr_t>(io.edgelist)));
    info("  edgemarkerlist = " +
         str(reinterpret_cast<std::uintptr_t>(io.edgemarkerlist)));
    info("  normlist = " + str(reinterpret_cast<std::uintptr_t>(io.normlist)));
    info("  numberofedges = " + str(io.numberofedges));
  }

  // Compute domain markers for subdomains
  static void compute_domain_markers(Mesh &mesh,
                                     const std::vector<Polygon> &subdomains)
  {
    info("Computing domain markers");
    Timer timer("compute_domain_markers");

    // build search tree for subdomains

    auto search_tree = BoundingBoxTree2D();
    std::vector<BoundingBox2D> bounding_boxes;
    for (const auto &subdomain : subdomains)
    {
      bounding_boxes.push_back(BoundingBox2D(subdomain));
    }
    search_tree.build(bounding_boxes);

    // Initialize domain markers and set all markers to -2 (ground)
    mesh.markers.resize(mesh.faces.size());
    std::fill(mesh.markers.begin(), mesh.markers.end(), -2);

    // Initialize markers for vertices belonging to a building
    std::vector<bool> is_building_vertex(mesh.vertices.size());
    std::fill(is_building_vertex.begin(), is_building_vertex.end(), false);

    // Iterate over cells to mark buildings
    if (subdomains.size() > 0)
    {
      for (size_t i = 0; i < mesh.faces.size(); i++)
      {
        // find building containg midpoint of cell (if any)
        const Vector3D c_3d = mesh.mid_point(i);
        const Vector2D c_2d(c_3d.x, c_3d.y);
        std::vector<size_t> indices = search_tree.find(Vector2D(c_2d));

        if (indices.size() > 0)
        {
          for (const auto &index : indices)
          {
            if (Geometry::polygon_contains_2d(subdomains[index], c_2d))
            {
              mesh.markers[i] = index;
              const Simplex2D &T = mesh.faces[i];
              // Mark all cell vertices as belonging to a building
              is_building_vertex[T.v0] = true;
              is_building_vertex[T.v1] = true;
              is_building_vertex[T.v2] = true;

              // // Check if individual vertices are inside a building
              // // (not only midpoint). Necessary for when building
              // // visualization meshes that are not boundary-fitted.
              // if (search_tree.find(mesh.vertices[T.v0]).size() == 0)
              //   is_building_vertex[T.v0] = false;
              // if (search_tree.find(mesh.vertices[T.v1]).size() == 0)
              //   is_building_vertex[T.v1] = false;
              // if (search_tree.find(mesh.vertices[T.v2]).size() == 0)
              //   is_building_vertex[T.v2] = false;

              break;
            }
          }
        }
      }

      // Iterate over cells to mark building halos
      for (size_t i = 0; i < mesh.faces.size(); i++)
      {
        // Check if any of the cell vertices belongs to a building
        const Simplex2D &T = mesh.faces[i];
        const bool touches_building =
            (is_building_vertex[T.v0] || is_building_vertex[T.v1] ||
             is_building_vertex[T.v2]);

        // Mark as halo (-1) if the cell touches a building but is not
        // itself inside footprint (not marked in the previous step)
        if (touches_building && mesh.markers[i] == -2)
          mesh.markers[i] = -1;
      }
    }
  }

  // Set x = min(x, y)
  static void set_min(double &x, double y)
  {
    if (y < x)
      x = y;
  }
};

} // namespace DTCC_BUILDER

#endif
