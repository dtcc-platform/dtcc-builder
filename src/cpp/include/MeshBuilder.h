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
  // Build mesh for city, returning a list of meshes.
  //
  // The first mesh is the ground (height map) and the remaining surfaces are
  // the extruded building footprints. Note that meshes are non-conforming; the
  // ground and building meshes are non-matching and the building meshes may
  // contain hanging nodes.
  static std::vector<Mesh>
  build_mesh(const City &city, const GridField &dtm, double resolution)
  {
    info("Building city mesh...");
    Timer timer("build_mesh");

    // Create empty subdomains for Triangle mesh building
    std::vector<std::vector<Point2D>> subDomains;

    // Get bounding box
    const BoundingBox2D &bbox = dtm.grid.BoundingBox;

    // Build boundary
    std::vector<Point2D> boundary;
    boundary.push_back(bbox.P);
    boundary.push_back(Point2D(bbox.Q.x, bbox.P.y));
    boundary.push_back(bbox.Q);
    boundary.push_back(Point2D(bbox.P.x, bbox.Q.y));

    // Build ground mesh
    Mesh ground_mesh;
    CallTriangle(ground_mesh, boundary, subDomains, resolution, false);

    // Compute domain markers
    ComputeDomainMarkers(ground_mesh, city);

    // Displace ground surface. Fill all points with maximum height. This is
    // used to always choose the smallest height for each point since each point
    // may be visited multiple times.
    const double zMax = dtm.Max();
    for (size_t i = 0; i < ground_mesh.Vertices.size(); i++)
      ground_mesh.Vertices[i].z = zMax;

    // If ground is not float, iterate over the triangles
    for (size_t i = 0; i < ground_mesh.Faces.size(); i++)
    {
      // Get cell marker
      const int cellMarker = ground_mesh.Markers[i];

      // Get triangle
      const Simplex2D &T = ground_mesh.Faces[i];

      // Check cell marker
      if (cellMarker != -2) // not ground
      {
        // Compute minimum height of vertices
        double zMin = std::numeric_limits<double>::max();
        zMin = std::min(zMin, dtm(ground_mesh.Vertices[T.v0]));
        zMin = std::min(zMin, dtm(ground_mesh.Vertices[T.v1]));
        zMin = std::min(zMin, dtm(ground_mesh.Vertices[T.v2]));

        // Set minimum height for all vertices
        setMin(ground_mesh.Vertices[T.v0].z, zMin);
        setMin(ground_mesh.Vertices[T.v1].z, zMin);
        setMin(ground_mesh.Vertices[T.v2].z, zMin);
      }
      else
      {
        // Sample height map at vertex position for all vertices
        setMin(ground_mesh.Vertices[T.v0].z, dtm(ground_mesh.Vertices[T.v0]));
        setMin(ground_mesh.Vertices[T.v1].z, dtm(ground_mesh.Vertices[T.v1]));
        setMin(ground_mesh.Vertices[T.v2].z, dtm(ground_mesh.Vertices[T.v2]));
      }
    }

    // Add ground mesh
    std::vector<Mesh> meshes;
    meshes.push_back(ground_mesh);

    // Get ground height (minimum)
    const double groundHeight = dtm.Min();

    // Iterate over buildings to build surfaces
    for (auto const &building : city.Buildings)
    {
      auto building_mesh = extrude_footprint(
          building.Footprint, resolution, groundHeight, building.MaxHeight());
      // Add surface
      meshes.push_back(building_mesh);
    }

    return meshes;
  }

  // Extrude Polygon to create a Mesh
  //
  static Mesh extrude_footprint(const Polygon &footprint,
                                double resolution,
                                double ground_height,
                                double height)
  {
    // FIXME: Consider making flipping triangles upside-down here
    // so that the normal points downwards rather than upwards.

    // Build 2D mesh of building footprint
    Mesh _mesh;
    // Create empty subdomains for Triangle mesh building
    // TODO: handle polygon with holes
    std::vector<std::vector<Point2D>> subDomains;

    CallTriangle(_mesh, footprint.Vertices, subDomains, resolution, false);
    // set ground height
    for (auto &v : _mesh.Vertices)
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
    const size_t numMeshPoints = _mesh.Vertices.size();
    const size_t numBoundaryPoints = footprint.Vertices.size();
    extrude_mesh.Vertices.resize(numMeshPoints + numBoundaryPoints);

    // Set total number of triangles
    const size_t numMeshTriangles = _mesh.Faces.size();
    const size_t numBoundaryTriangles = 2 * numBoundaryPoints;
    extrude_mesh.Faces.resize(numMeshTriangles + numBoundaryTriangles);

    // Add points at top
    for (size_t i = 0; i < numMeshPoints; i++)
    {
      const Point3D &p2D = _mesh.Vertices[i];
      const Vector3D p3D(p2D.x, p2D.y, height);
      extrude_mesh.Vertices[i] = p3D;
    }

    // Add points at bottom
    for (size_t i = 0; i < numBoundaryPoints; i++)
    {
      const Point3D &p2D = _mesh.Vertices[i];
      const Vector3D p3D(p2D.x, p2D.y, ground_height);
      extrude_mesh.Vertices[numMeshPoints + i] = p3D;
    }

    // Add triangles on top
    for (size_t i = 0; i < numMeshTriangles; i++)
      extrude_mesh.Faces[i] = _mesh.Faces[i];

    // Add triangles on boundary
    for (size_t i = 0; i < numBoundaryPoints; i++)
    {
      const size_t v0 = i;
      const size_t v1 = (i + 1) % numBoundaryPoints;
      const size_t v2 = v0 + numMeshPoints;
      const size_t v3 = v1 + numMeshPoints;
      Simplex2D t0(v0, v2, v1); // Outward-pointing normal
      Simplex2D t1(v1, v2, v3); // Outward-pointing normal
      extrude_mesh.Faces[numMeshTriangles + 2 * i] = t0;
      extrude_mesh.Faces[numMeshTriangles + 2 * i + 1] = t1;
    }

    return extrude_mesh;
  }

  // Build ground mesh for city.
  //
  // The mesh is a triangular mesh of the rectangular region
  // defined by (xmin, xmax) x (ymin, ymax). The edges of the mesh respect
  // the boundaries of the buildings.
  //
  // Markers:
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

    // Print some stats
    const BoundingBox2D boundingBox(Point2D(xmin, ymin), Point2D(xmax, ymax));
    const size_t nx = (boundingBox.Q.x - boundingBox.P.x) / resolution;
    const size_t ny = (boundingBox.Q.y - boundingBox.P.y) / resolution;
    const size_t n = nx * ny;
    info("Domain bounding box is " + str(boundingBox));
    info("Mesh resolution is " + str(resolution));
    info("Estimated number of triangles is " + str(n));
    info("Number of subdomains (buildings) is " + str(city.Buildings.size()));

    // Extract subdomains (building footprints)
    std::vector<std::vector<Point2D>> subDomains;
    for (auto const &building : city.Buildings)
      subDomains.push_back(building.Footprint.Vertices);

    // Build boundary
    std::vector<Point2D> boundary{};
    boundary.push_back(boundingBox.P);
    boundary.push_back(Point2D(boundingBox.Q.x, boundingBox.P.y));
    boundary.push_back(boundingBox.Q);
    boundary.push_back(Point2D(boundingBox.P.x, boundingBox.Q.y));

    // Build 2D mesh
    Mesh mesh;
    CallTriangle(mesh, boundary, subDomains, resolution, true);

    // Mark subdomains
    ComputeDomainMarkers(mesh, city);

    return mesh;
  }

  // Layer ground mesh to create a volume mesh.
  //
  // The volume mesh is a tetrahedral mesh constructed
  // extruding the 2D mesh in the vertical (z) direction.
  //
  // Markers:
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
                                      double domainHeight,
                                      double meshResolution)
  {
    Timer timer("build_volume_mesh");

    // Compute number of layers
    const size_t numLayers = int(std::ceil(domainHeight / meshResolution));
    const double dz = domainHeight / double(numLayers);
    const size_t layerSize = ground_mesh.Vertices.size();

    info("Building volume mesh with " + str(numLayers) + " layers...");

    // Initialize volume mesh
    VolumeMesh volume_mesh;
    volume_mesh.Vertices.resize((numLayers + 1) * ground_mesh.Vertices.size());
    volume_mesh.Cells.resize(numLayers * 3 * ground_mesh.Faces.size());
    volume_mesh.Markers.resize(volume_mesh.Cells.size());
    volume_mesh.num_layers = numLayers;

    // Add vertices
    {
      size_t k = 0;
      for (size_t layer = 0; layer <= numLayers; layer++)
      {
        // Compute height of layer
        const double z = layer * dz;

        // Iterate over vertices in layer
        for (const auto &p2D : ground_mesh.Vertices)
          volume_mesh.Vertices[k++] = Point3D(p2D.x, p2D.y, z);
      }
    }

    // Add cells
    {
      size_t k = 0;
      size_t offset = 0;
      for (size_t layer = 0; layer < numLayers; layer++)
      {
        // Iterate over triangles in layer
        for (const auto &T : ground_mesh.Faces)
        {
          // Get sorted vertex indices for bottom layer
          const size_t u0 = T.v0 + offset;
          const size_t u1 = T.v1 + offset;
          const size_t u2 = T.v2 + offset;

          // Get sorted vertices for top layer
          const size_t v0 = u0 + layerSize;
          const size_t v1 = u1 + layerSize;
          const size_t v2 = u2 + layerSize;

          // Create three tetrahedra by connecting the first vertex
          // of each edge in the bottom layer with the second
          // vertex of the corresponding edge in the top layer.
          volume_mesh.Cells[k++] = Simplex3D(u0, u1, u2, v2);
          volume_mesh.Cells[k++] = Simplex3D(u0, v1, u1, v2);
          volume_mesh.Cells[k++] = Simplex3D(u0, v0, v1, v2);
        }

        // Add to offset
        offset += layerSize;
      }
    }

    // Add domain markers
    {
      size_t k = 0;
      for (size_t layer = 0; layer < numLayers; layer++)
      {
        for (const auto &marker : ground_mesh.Markers)
        {
          int m = 0;

          // Top layer marked as -3
          if (layer == numLayers - 1)
            m = -3;

          // Halo and ground only marked for bottom layer
          else if (marker == -1 || marker == -2)
          {
            if (layer == 0)
              m = marker;
            else
              m = -4;
          }

          // Buildings marked for all layers (except top layer).
          // Later adjusted to -4 above buildings in trim_volume_mesh.
          else
          {
            m = marker;
          }

          // Set markers
          volume_mesh.Markers[k++] = m;
          volume_mesh.Markers[k++] = m;
          volume_mesh.Markers[k++] = m;
        }
      }
    }

    return volume_mesh;
  }

  // Trim volume mesh by removing cells inside buildings.
  //
  // Markers:
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
  // - Cells in top layer marked as top (-3)
  // - All other cells (in between) marked as other (-4)
  static VolumeMesh trim_volume_mesh(const VolumeMesh &volume_mesh,
                                     const Mesh &mesh,
                                     const City &city)
  {
    info("Trimming volume mesh...");
    Timer timer("trim_volume_mesh");

    // Get sizes
    const size_t numBuildings = city.Buildings.size();
    const size_t numCells2D = mesh.Faces.size();
    const size_t numCells3D = volume_mesh.Cells.size();
    const size_t layerSize = 3 * mesh.Faces.size();

    // Phase 1: Determine which cells should be trimmed
    // ------------------------------------------------

    // Build map from buildings to cells in 2D mesh
    std::vector<std::vector<size_t>> buildingCells2D(numBuildings);
    for (size_t cellIndex2D = 0; cellIndex2D < numCells2D; cellIndex2D++)
    {
      const int buildingIndex = mesh.Markers[cellIndex2D];
      if (buildingIndex >= 0)
        buildingCells2D[buildingIndex].push_back(cellIndex2D);
    }

    // Create markers for cells to be trimmed (keep by default)
    std::vector<bool> trimCell(numCells3D);
    std::fill(trimCell.begin(), trimCell.end(), false);

    // Keep track of first layer for each building
    std::vector<size_t> firstLayer(numBuildings);
    std::fill(firstLayer.begin(), firstLayer.end(), 0);

    // Iterate over buildings
    for (size_t buildingIndex = 0; buildingIndex < numBuildings;
         buildingIndex++)
    {
      // Iterate over layers
      for (size_t layer = 0; layer < volume_mesh.num_layers; layer++)
      {
        // Build list of 3D cells for building in current layer
        std::vector<size_t> cells3D;
        for (const auto &cellIndex2D : buildingCells2D[buildingIndex])
        {
          for (size_t j = 0; j < 3; j++)
            cells3D.push_back(Index3D(layer, layerSize, cellIndex2D, j));
        }

        // Trim layer if any cell midpoint is below building height
        bool trimLayer = false;
        for (const auto &cellIndex3D : cells3D)
        {
          const double z = volume_mesh.MidPoint(cellIndex3D).z;
          const double h = city.Buildings[buildingIndex].MaxHeight();
          if (z < h)
          {
            trimLayer = true;
            break;
          }
        }

        // Check if layer should be trimmed
        if (trimLayer)
        {
          // Mark cells for trimming
          for (const auto &cellIndex3D : cells3D)
            trimCell[cellIndex3D] = true;
        }
        else
        {
          // If layer should be kept, no need to check more layers
          firstLayer[buildingIndex] = layer;
          break;
        }
      }
    }

    // Phase 2: Adjust markers
    // -----------------------

    // Create copy of markeres
    std::vector<int> markers{volume_mesh.Markers};

    // Mark cells between bottom and top layer as -4
    for (size_t layer = 1; layer < volume_mesh.num_layers - 1; layer++)
    {
      for (size_t cellIndex2D = 0; cellIndex2D < numCells2D; cellIndex2D++)
        for (size_t j = 0; j < 3; j++)
          markers[Index3D(layer, layerSize, cellIndex2D, j)] = -4;
    }

    // Mark cells in top layer as -3
    for (size_t cellIndex2D = 0; cellIndex2D < numCells2D; cellIndex2D++)
    {
      for (size_t j = 0; j < 3; j++)
        markers[Index3D(volume_mesh.num_layers - 1, layerSize, cellIndex2D,
                        j)] = -3;
    }

    // Mark cells in first layer above each building:
    //
    // 0, 1, 2, ... if building is not covered by bottom layer (normal case)
    // -1           if building is covered by bottom layer (modify to halo)
    for (size_t buildingIndex = 0; buildingIndex < numBuildings;
         buildingIndex++)
    {
      const size_t layer = firstLayer[buildingIndex];
      size_t marker = buildingIndex;
      if (layer == 0)
      {
        warning("Building " + str(buildingIndex) +
                " is covered by bottom layer");
        marker = -1;
      }
      for (const auto &cellIndex2D : buildingCells2D[buildingIndex])
      {
        for (size_t j = 0; j < 3; j++)
          markers[Index3D(layer, layerSize, cellIndex2D, j)] = marker;
      }
    }

    // Phase 3: Extract new mesh for all cells that should be kept
    // -----------------------------------------------------------

    // Renumber vertices and cells
    std::unordered_map<size_t, size_t> vertexMap;
    std::unordered_map<size_t, size_t> cellMap;
    size_t k = 0;
    size_t l = 0;
    for (size_t cellIndex3D = 0; cellIndex3D < volume_mesh.Cells.size();
         cellIndex3D++)
    {
      if (!trimCell[cellIndex3D])
      {
        // Get cell
        const Simplex3D &T = volume_mesh.Cells[cellIndex3D];

        // Renumbers vertices
        if (vertexMap.find(T.v0) == vertexMap.end())
          vertexMap[T.v0] = k++;
        if (vertexMap.find(T.v1) == vertexMap.end())
          vertexMap[T.v1] = k++;
        if (vertexMap.find(T.v2) == vertexMap.end())
          vertexMap[T.v2] = k++;
        if (vertexMap.find(T.v3) == vertexMap.end())
          vertexMap[T.v3] = k++;

        // Renumber cells
        cellMap[cellIndex3D] = l++;
      }
    }

    // Initialize new mesh data
    const size_t numVertices = vertexMap.size();
    const size_t numCells = cellMap.size();
    std::vector<Point3D> _vertices(numVertices);
    std::vector<Simplex3D> _cells(numCells);
    std::vector<int> _markers(numCells);

    // Set new mesh data
    for (const auto v : vertexMap)
      _vertices[v.second] = volume_mesh.Vertices[v.first];
    for (const auto c : cellMap)
    {
      _cells[c.second].v0 = vertexMap[volume_mesh.Cells[c.first].v0];
      _cells[c.second].v1 = vertexMap[volume_mesh.Cells[c.first].v1];
      _cells[c.second].v2 = vertexMap[volume_mesh.Cells[c.first].v2];
      _cells[c.second].v3 = vertexMap[volume_mesh.Cells[c.first].v3];
      _markers[c.second] = markers[c.first];
    }

    // Create new mesh and assign data
    VolumeMesh _volume_mesh;
    _volume_mesh.Vertices = _vertices;
    _volume_mesh.Cells = _cells;
    _volume_mesh.Markers = _markers;

    return _volume_mesh;
  }

private:
  // Map from 2D cell index to 3D cell indices
  static size_t
  Index3D(size_t layer, size_t layerSize, size_t cellIndex2D, size_t j)
  {
    return layer * layerSize + 3 * cellIndex2D + j;
  }

  // Call Triangle to compute 2D mesh
  static void CallTriangle(Mesh &mesh,
                           const std::vector<Point2D> &boundary,
                           const std::vector<std::vector<Point2D>> &subDomains,
                           double h,
                           bool sortTriangles)
  {
    Timer timer("CallTriangle");

    // Set area constraint to control mesh size
    const double maxArea = 0.5 * h * h;

    // Set input switches for Triangle
    char triswitches[64];
    // sprintf(triswitches, "zQpq25a%.16f", maxArea);
    snprintf(triswitches, sizeof(triswitches), "zQpq25a%.16f", maxArea);

    // z = use zero-based numbering
    // p = use polygon input (segments)
    // q = control mesh quality
    //
    // Note that the minimum angle (here 25) should be
    // as large as possible for high quality meshes but
    // it should be less than 28.6 degrees to guarantee
    // that Triangle terminates. Default is 20 degrees.

    // Create input data structure for Triangle
    struct triangulateio in = CreateTriangleIO();

    // Set number of points
    size_t numPoints = boundary.size();
    for (auto const &innerPolygon : subDomains)
      numPoints += innerPolygon.size();
    in.numberofpoints = numPoints;

    // Set points
    in.pointlist = new double[2 * numPoints];
    {
      size_t k = 0;
      for (auto const &p : boundary)
      {
        in.pointlist[k++] = p.x;
        in.pointlist[k++] = p.y;
      }
      for (auto const &innerPolygon : subDomains)
      {
        for (auto const &p : innerPolygon)
        {
          in.pointlist[k++] = p.x;
          in.pointlist[k++] = p.y;
        }
      }
    }

    // Set number of segments
    const size_t numSegments = numPoints;
    in.numberofsegments = numSegments;

    // Set segments
    in.segmentlist = new int[2 * numSegments];
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
      for (size_t i = 0; i < subDomains.size(); i++)
      {
        for (size_t j = 0; j < subDomains[i].size(); j++)
        {
          const size_t j0 = j;
          const size_t j1 = (j + 1) % subDomains[i].size();
          in.segmentlist[k++] = n + j0;
          in.segmentlist[k++] = n + j1;
        }
        n += subDomains[i].size();
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
    Point2D c;
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
    struct triangulateio out = CreateTriangleIO();
    struct triangulateio vorout = CreateTriangleIO();

    // Call Triangle
    triangulate(triswitches, &in, &out, &vorout);

    // Uncomment for debugging
    // PrintTriangleIO(out);
    // PrintTriangleIO(vorout);

    // Extract points
    mesh.Vertices.reserve(out.numberofpoints);
    for (int i = 0; i < out.numberofpoints; i++)
    {
      Point3D p(out.pointlist[2 * i], out.pointlist[2 * i + 1], 0.0);
      mesh.Vertices.push_back(p);
    }

    // Extract triangles
    mesh.Faces.reserve(out.numberoftriangles);
    for (int i = 0; i < out.numberoftriangles; i++)
    {
      // Note the importance of creating a sorted simplex here!
      Simplex2D t(out.trianglelist[3 * i], out.trianglelist[3 * i + 1],
                  out.trianglelist[3 * i + 2], sortTriangles);
      mesh.Faces.push_back(t);
    }

    // Free memory
    // trifree(&out); // causes segfault
    delete[] in.pointlist;
    delete[] in.segmentlist;
    delete[] in.holelist;
  }

  // Create and reset Triangle I/O data structure
  static struct triangulateio CreateTriangleIO()
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

  // Print triangle I/O data
  static void PrintTriangleIO(const struct triangulateio &io)
  {
    info("Triangle I/O data: ");
    info("  pointlist = " +
         str(reinterpret_cast<std::uintptr_t>(io.pointlist)));
    info("  pointmarkerlist = " +
         str(reinterpret_cast<std::uintptr_t>(io.pointmarkerlist)));
    if (io.pointmarkerlist)
    {
      std::stringstream stringBuilder{};
      stringBuilder << "   ";
      for (int i = 0; i < io.numberofpoints; i++)
        stringBuilder << " " << io.pointmarkerlist[i];
      stringBuilder << std::endl;
      info(stringBuilder.str());
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
      std::stringstream stringBuilder{};
      stringBuilder << "   ";
      for (int i = 0; i < io.numberofsegments; i++)
        stringBuilder << " " << io.segmentmarkerlist[i];
      stringBuilder << std::endl;
      info(stringBuilder.str());
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
  static void ComputeDomainMarkers(Mesh &mesh, const City &city)
  {
    info("Computing domain markers");
    Timer timer("ComputeDomainMarkers");

    // Build search tree for city
    city.BuildSearchTree();

    // Initialize domain markers and set all markers to -2 (ground)
    mesh.Markers.resize(mesh.Faces.size());
    std::fill(mesh.Markers.begin(), mesh.Markers.end(), -2);

    // Initialize markers for vertices belonging to a building
    std::vector<bool> isBuildingVertex(mesh.Vertices.size());
    std::fill(isBuildingVertex.begin(), isBuildingVertex.end(), false);

    // Iterate over cells to mark buildings
    for (size_t i = 0; i < mesh.Faces.size(); i++)
    {
      // Find building containg midpoint of cell (if any)
      const Point3D c3d = mesh.MidPoint(i);
      const Point2D c2d(c3d.x, c3d.y);
      const int marker = city.FindBuilding(Vector2D(c2d));

      // Get triangle
      const Simplex2D &T = mesh.Faces[i];

      // Check if we are inside a building
      if (marker >= 0)
      {
        // Set domain marker to building number
        mesh.Markers[i] = marker;

        // Mark all cell vertices as belonging to a building
        isBuildingVertex[T.v0] = true;
        isBuildingVertex[T.v1] = true;
        isBuildingVertex[T.v2] = true;
      }

      // Check if individual vertices are inside a building
      // (not only midpoint). Necessary for when building
      // visualization meshes that are not boundary-fitted.
      if (city.FindBuilding(Vector3D(mesh.Vertices[T.v0])) >= 0)
        isBuildingVertex[T.v0] = true;
      if (city.FindBuilding(Vector3D(mesh.Vertices[T.v1])) >= 0)
        isBuildingVertex[T.v1] = true;
      if (city.FindBuilding(Vector3D(mesh.Vertices[T.v2])) >= 0)
        isBuildingVertex[T.v2] = true;
    }

    // Iterate over cells to mark building halos
    for (size_t i = 0; i < mesh.Faces.size(); i++)
    {
      // Check if any of the cell vertices belongs to a building
      const Simplex2D &T = mesh.Faces[i];
      const bool touchesBuilding =
          (isBuildingVertex[T.v0] || isBuildingVertex[T.v1] ||
           isBuildingVertex[T.v2]);

      // Mark as halo (-1) if the cell touches a building but is not
      // itself inside footprint (not marked in the previous step)
      if (touchesBuilding && mesh.Markers[i] == -2)
        mesh.Markers[i] = -1;
    }
  }

  // Set x = min(x, y)
  static void setMin(double &x, double y)
  {
    if (y < x)
      x = y;
  }
};

} // namespace DTCC_BUILDER

#endif
