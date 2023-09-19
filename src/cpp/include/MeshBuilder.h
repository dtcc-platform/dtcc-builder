// Copyright (C) 2018 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_BUILDER_H
#define DTCC_MESH_BUILDER_H

#include <cmath>
#include <iostream>
#include <stack>
#include <tuple>
#include <vector>

#include "Geometry.h"
#include "Logging.h"
#include "Timer.h"
#include "model/City.h"
#include "model/GridField.h"
#include "model/Mesh.h"
#include "model/Vector.h"

extern "C"
{
#include <triangle.h>
}

namespace DTCC_BUILDER
{

class MeshBuilder
{
public:
  // build mesh for city, returning a list of meshes.
  //
  // The first mesh is the ground (height map) and the remaining surfaces are
  // the extruded building footprints. Note that meshes are non-conforming; the
  // ground and building meshes are non-matching and the building meshes may
  // contain hanging nodes.
  static std::vector<Mesh> build_mesh(const City &city,
                                      const GridField &dtm,
                                      double resolution,
                                      bool ground_only = false)
  {
    info("Building city mesh...");
    Timer timer("build_mesh");

    // Create empty subdomains for Triangle mesh building
    std::vector<std::vector<Vector2D>> sub_domains;

    // Get bounding box
    const BoundingBox2D &bbox = dtm.grid.bounding_box;

    // build boundary
    std::vector<Vector2D> boundary;
    boundary.push_back(bbox.P);
    boundary.push_back(Vector2D(bbox.Q.x, bbox.P.y));
    boundary.push_back(bbox.Q);
    boundary.push_back(Vector2D(bbox.P.x, bbox.Q.y));

    // build ground mesh
    Mesh ground_mesh;
    call_triangle(ground_mesh, boundary, sub_domains, resolution, false);

    // Compute domain markers
    compute_domain_markers(ground_mesh, city);

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

    // Add ground mesh
    std::vector<Mesh> meshes;
    meshes.push_back(ground_mesh);

    // Get ground height (minimum)
    const double ground_height = dtm.min();

    if (!ground_only)
    {
      // Iterate over buildings to build surfaces
      for (auto const &building : city.buildings)
      {
        auto building_mesh =
            extrude_footprint(building.footprint, resolution, ground_height,
                              building.max_height());
        // Add surface
        meshes.push_back(building_mesh);
      }
    }
    return meshes;
  }

  // static Mesh
  // build_ground_mesh(const City &city, const GridField &dtm, double
  // resolution)
  // {
  // }

  // Extrude Polygon to create a Mesh
  //
  static Mesh extrude_footprint(const Polygon &footprint,
                                double resolution,
                                double ground_height,
                                double height,
                                bool cap_base = false)
  {
    // FIXME: Consider making flipping triangles upside-down here
    // so that the normal points downwards rather than upwards.

    // build 2D mesh of building footprint
    Mesh _mesh;
    // Create empty subdomains for Triangle mesh building
    // TODO: handle polygon with holes
    std::vector<std::vector<Vector2D>> sub_domains;

    call_triangle(_mesh, footprint.vertices, sub_domains, resolution, false);
    // set ground height
    for (auto &v : _mesh.vertices)
      v.z = ground_height;

    // Create empty building surface
    Mesh extrude_mesh;

    // Note: The 2D mesh contains all the input boundary points with
    // the same numbers as in the footprint polygon, but may also
    // contain new points (Steiner points) added during mesh
    // generation. We add the top points (including any Steiner
    // points) first, then the points at the bottom (the footprint).

    // Get absolute height of building
    // const double buildingHeight = height;

    // Set total number of points
    const size_t num_mesh_points = _mesh.vertices.size();
    const size_t num_boundary_points = footprint.vertices.size();
    size_t total_points = num_mesh_points + num_boundary_points;
    if (cap_base)
      total_points += num_mesh_points; // add points for base cap
    extrude_mesh.vertices.resize(total_points);

    // Set total number of triangles
    const size_t num_mesh_triangles = _mesh.faces.size();
    const size_t num_boundary_triangles = 2 * num_boundary_points;
    size_t total_faces = num_mesh_triangles + num_boundary_triangles;
    if (cap_base)
      total_faces += num_mesh_triangles; // add triangles for base cap
    extrude_mesh.faces.resize(total_faces);

    // Add points at top
    for (size_t i = 0; i < num_mesh_points; i++)
    {
      const Vector3D &p_2d = _mesh.vertices[i];
      const Vector3D p_3d(p_2d.x, p_2d.y, height);
      extrude_mesh.vertices[i] = p_3d;
    }

    // Add points at bottom
    for (size_t i = 0; i < num_boundary_points; i++)
    {
      const Vector3D &p_2d = _mesh.vertices[i];
      const Vector3D p_3d(p_2d.x, p_2d.y, ground_height);
      extrude_mesh.vertices[num_mesh_points + i] = p_3d;
    }

    // Add triangles on top
    for (size_t i = 0; i < num_mesh_triangles; i++)
      extrude_mesh.faces[i] = _mesh.faces[i];

    // Add triangles on boundary
    for (size_t i = 0; i < num_boundary_points; i++)
    {
      const size_t v0 = i;
      const size_t v1 = (i + 1) % num_boundary_points;
      const size_t v2 = v0 + num_mesh_points;
      const size_t v3 = v1 + num_mesh_points;
      Simplex2D t0(v0, v2, v1); // Outward-pointing normal
      Simplex2D t1(v1, v2, v3); // Outward-pointing normal
      extrude_mesh.faces[num_mesh_triangles + 2 * i] = t0;
      extrude_mesh.faces[num_mesh_triangles + 2 * i + 1] = t1;
    }

    if (cap_base)
    {
      // Add points for base
      const size_t vertex_offset = num_mesh_points + num_boundary_points;
      for (size_t i = 0; i < num_mesh_points; i++)
      {
        const Vector3D &p2D = _mesh.vertices[i];
        const Vector3D p3D(p2D.x, p2D.y, ground_height);
        extrude_mesh.vertices[vertex_offset + i] = p3D;
      }
      // Add triangles on top
      const size_t face_offset = num_mesh_triangles + num_boundary_triangles;
      for (size_t i = 0; i < num_mesh_triangles; i++)
      {
        auto face = _mesh.faces[i];
        face.v0 += vertex_offset;
        face.v1 += vertex_offset;
        face.v2 += vertex_offset;
        std::swap(face.v1, face.v2); // flip triangle
        extrude_mesh.faces[face_offset + i] = face;
      }
    }

    return extrude_mesh;
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
  static Mesh build_ground_mesh(const City &city,
                                double xmin,
                                double ymin,
                                double xmax,
                                double ymax,
                                double resolution)
  {
    info("Building ground mesh for city...");
    Timer timer("build_ground_mesh");

    // print some stats
    const BoundingBox2D bounding_box(Vector2D(xmin, ymin),
                                     Vector2D(xmax, ymax));
    const size_t nx = (bounding_box.Q.x - bounding_box.P.x) / resolution;
    const size_t ny = (bounding_box.Q.y - bounding_box.P.y) / resolution;
    const size_t n = nx * ny;
    info("Domain bounding box is " + str(bounding_box));
    info("Mesh resolution is " + str(resolution));
    info("Estimated number of triangles is " + str(n));
    info("Number of subdomains (buildings) is " + str(city.buildings.size()));

    // Extract subdomains (building footprints)
    std::vector<std::vector<Vector2D>> sub_domains;
    for (auto const &building : city.buildings)
      sub_domains.push_back(building.footprint.vertices);

    // build boundary
    std::vector<Vector2D> boundary{};
    boundary.push_back(bounding_box.P);
    boundary.push_back(Vector2D(bounding_box.Q.x, bounding_box.P.y));
    boundary.push_back(bounding_box.Q);
    boundary.push_back(Vector2D(bounding_box.P.x, bounding_box.Q.y));

    // build 2D mesh
    Mesh mesh;
    call_triangle(mesh, boundary, sub_domains, resolution, true);

    // Mark subdomains
    compute_domain_markers(mesh, city);

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
                                      double mesh_resolution)
  {
    Timer timer("build_volume_mesh");

    // Compute number of layers
    const size_t num_layers = int(std::ceil(domain_height / mesh_resolution));
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
                                     const City &city)
  {
    info("Trimming volume mesh...");
    Timer timer("trim_volume_mesh");

    // Get sizes
    const size_t num_buildings = city.buildings.size();
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
          const double h = city.buildings[building_index].max_height();
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

private:
  // Map from 2D cell index to 3D cell indices
  static size_t
  index_3d(size_t layer, size_t layer_size, size_t cell_index_2d, size_t j)
  {
    return layer * layer_size + 3 * cell_index_2d + j;
  }

  // Call Triangle to compute 2D mesh
  static void
  call_triangle(Mesh &mesh,
                const std::vector<Vector2D> &boundary,
                const std::vector<std::vector<Vector2D>> &sub_domains,
                double h,
                bool sort_triangles)
  {
    Timer timer("call_triangle");

    // Set area constraint to control mesh size
    const double max_area = 0.5 * h * h;

    // Set input switches for Triangle
    char triswitches[64];
    // sprintf(triswitches, "zQpq25a%.16f", max_area);
    snprintf(triswitches, sizeof(triswitches), "zQpq25a%.16f", max_area);

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
      num_points += innerPolygon.size();
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
    triangulate(triswitches, &in, &out, &vorout);

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
  static void compute_domain_markers(Mesh &mesh, const City &city)
  {
    info("Computing domain markers");
    Timer timer("compute_domain_markers");

    // build search tree for city
    city.build_search_tree();

    // Initialize domain markers and set all markers to -2 (ground)
    mesh.markers.resize(mesh.faces.size());
    std::fill(mesh.markers.begin(), mesh.markers.end(), -2);

    // Initialize markers for vertices belonging to a building
    std::vector<bool> is_building_vertex(mesh.vertices.size());
    std::fill(is_building_vertex.begin(), is_building_vertex.end(), false);

    // Iterate over cells to mark buildings
    for (size_t i = 0; i < mesh.faces.size(); i++)
    {
      // find building containg midpoint of cell (if any)
      const Vector3D c_3d = mesh.mid_point(i);
      const Vector2D c_2d(c_3d.x, c_3d.y);
      const int marker = city.find_building(Vector2D(c_2d));

      // Get triangle
      const Simplex2D &T = mesh.faces[i];

      // Check if we are inside a building
      if (marker >= 0)
      {
        // Set domain marker to building number
        mesh.markers[i] = marker;

        // Mark all cell vertices as belonging to a building
        is_building_vertex[T.v0] = true;
        is_building_vertex[T.v1] = true;
        is_building_vertex[T.v2] = true;
      }

      // Check if individual vertices are inside a building
      // (not only midpoint). Necessary for when building
      // visualization meshes that are not boundary-fitted.
      if (city.find_building(Vector3D(mesh.vertices[T.v0])) >= 0)
        is_building_vertex[T.v0] = true;
      if (city.find_building(Vector3D(mesh.vertices[T.v1])) >= 0)
        is_building_vertex[T.v1] = true;
      if (city.find_building(Vector3D(mesh.vertices[T.v2])) >= 0)
        is_building_vertex[T.v2] = true;
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

  // Set x = min(x, y)
  static void set_min(double &x, double y)
  {
    if (y < x)
      x = y;
  }
};

} // namespace DTCC_BUILDER

#endif
