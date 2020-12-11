// Copyright (C) 2018 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_GENERATOR_H
#define DTCC_MESH_GENERATOR_H

#include <cmath>
#include <iostream>
#include <stack>
#include <tuple>
#include <vector>

#include "CSV.h"
#include "CityModel.h"
#include "Geometry.h"
#include "Mesh.h"
#include "Vector.h"
#include "Surface.h"
#include "GridField.h"
#include "Timer.h"
#include "Logging.h"

extern "C"
{
#include <triangle.h>
}

namespace DTCC
{

  class MeshGenerator
  {
  public:
    // Generate 2D mesh. The mesh is a triangular mesh of the rectangular region
    // defined by (xMin, xMax) x (yMin, yMax). The edges of the mesh respect the
    // boundaries of the buildings.
    //
    // The domain markers label the triangles inside building footprints
    // with the number of the building (0, 1, 2, ...). Triangles that are
    // neighbors of building triangles (just outside buildings) are marked
    // as -1 and the remaining triangles (the ground) are marked as -2.
    static void GenerateMesh2D(Mesh2D &mesh2D,
                               const CityModel &cityModel,
                               const BoundingBox2D &boundingBox,
                               double resolution)
    {
      Info("MeshGenerator: Generating 2D mesh...");
      Timer("GenerateMesh2D");

      // Extract subdomains (building footprints)
      std::vector<std::vector<Point2D>> subDomains;
      for (auto const &building : cityModel.Buildings)
        subDomains.push_back(building.Footprint.Vertices);

      // Generate boundary
      std::vector<Point2D> boundary{};
      boundary.push_back(boundingBox.P);
      boundary.push_back(Point2D(boundingBox.Q.x, boundingBox.P.y));
      boundary.push_back(boundingBox.Q);
      boundary.push_back(Point2D(boundingBox.P.x, boundingBox.Q.y));

      // Generate 2D mesh
      CallTriangle(mesh2D, boundary, subDomains, resolution);

      // Mark subdomains
      ComputeDomainMarkers(mesh2D, cityModel);
    }

    // Generate 3D mesh. The mesh is a tetrahedral mesh generated by
    // extruding the 2D mesh in the vertical (z) direction.
    //
    // The domain markers for the 3D mesh are propagated as is from
    // the domain markers of the 2D mesh, except for cells in the top
    // layer that are marked as -3.
    //
    // It is assumed that the 2D mesh is sorted, which is the case if the
    // mesh has been generated by calling GenerateMesh2D().
    //
    // @return Number of layers
    static size_t GenerateMesh3D(Mesh3D &mesh3D,
                                 const Mesh2D &mesh2D,
                                 double domainHeight,
                                 double meshResolution)
    {
      Info("MeshGenerator: Generating 3D mesh...");
      Timer("GenerateMesh3D");

      // Compute number of layers
      const size_t numLayers = int(std::ceil(domainHeight / meshResolution));
      const double dz = domainHeight / double(numLayers);
      const size_t layerSize = mesh2D.Vertices.size();

      Info("MeshGenerator: Generating 3D mesh with " + str(numLayers) +
           " layers...");

      // Resize 3D mesh
      mesh3D.Vertices.resize((numLayers + 1) * mesh2D.Vertices.size());
      mesh3D.Cells.resize(numLayers * 3 * mesh2D.Cells.size());
      mesh3D.Markers.resize(mesh3D.Cells.size());

      // Add vertices
      {
        size_t k = 0;
        for (size_t layer = 0; layer <= numLayers; layer++)
        {
          // Compute height of layer
          const double z = layer * dz;

          // Iterate over vertices in layer
          for (const auto &p2D : mesh2D.Vertices)
            mesh3D.Vertices[k++] = Point3D(p2D.x, p2D.y, z);
        }
      }

      // Add cells
      {
        size_t k = 0;
        size_t offset = 0;
        for (size_t layer = 0; layer < numLayers; layer++)
        {
          // Iterate over triangles in layer
          for (const auto &T : mesh2D.Cells)
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
            mesh3D.Cells[k++] = Simplex3D(u0, u1, u2, v2);
            mesh3D.Cells[k++] = Simplex3D(u0, v1, u1, v2);
            mesh3D.Cells[k++] = Simplex3D(u0, v0, v1, v2);
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
          for (const auto &marker : mesh2D.Markers)
          {
            const int m = (layer == numLayers - 1 ? -3 : marker);
            mesh3D.Markers[k++] = m;
            mesh3D.Markers[k++] = m;
            mesh3D.Markers[k++] = m;
          }
        }
      }

      return numLayers;
    }

    // Trim 3D mesh. The mesh is trimmed by removing tetrahedra inside
    // building shape. The domain markers are left untouched.
    static void TrimMesh3D(Mesh3D &mesh3D,
                           const Mesh2D &mesh2D,
                           const CityModel &cityModel,
                           size_t numLayers)
    {
      Info("MeshGenerator: Trimming 3D mesh...");
      Timer("TrimMesh3D");

      // Create markers for all cells that should be kept
      std::vector<bool> keepCell(mesh3D.Cells.size());
      std::fill(keepCell.begin(), keepCell.end(), true);

      // Create markers for all buildings that should be kept in current layer
      std::vector<bool> keepBuilding(cityModel.Buildings.size());

      // Phase 1: Mark which cells that should be kept. This requires some
      // care since we want to mark all cells in a layer that belong to
      // the same building in the same way.

      // Iterate over layers
      const size_t layerSize = 3 * mesh2D.Cells.size();
      for (size_t layer = 0; layer < numLayers; layer++)
      {
        // Keep all buildings by default in current layer
        std::fill(keepBuilding.begin(), keepBuilding.end(), true);

        // Iterate over triangles in 2D mesh
        for (size_t i = 0; i < mesh2D.Cells.size(); i++)
        {
          // Get marker (number of building)
          const int marker = mesh2D.Markers[i];

          // Check if we're inside a building
          if (marker >= 0)
          {
            // Iterate over building cells in current layer
            for (size_t j = 0; j < 3; j++)
            {
              // Get index of cell
              const size_t cellIndex = layer * layerSize + 3 * i + j;

              // Check if midpoint is below building height
              Point3D p = mesh3D.MidPoint(cellIndex);
              if (p.z < cityModel.Buildings[marker].MaxHeight())
              {
                keepBuilding[marker] = false;
                break;
              }
            }
          }
        }

        // Now that we know how to treat each building, iterate
        // over the cells again to mark them for removal

        // Iterate over triangles in 2D mesh
        for (size_t i = 0; i < mesh2D.Cells.size(); i++)
        {
          // Get marker (number of building)
          const int marker = mesh2D.Markers[i];

          // Check if we're inside a building
          if (marker >= 0)
          {
            // Iterate over building cells in current layer
            for (size_t j = 0; j < 3; j++)
            {
              // Get index of cell
              const size_t cellIndex = layer * layerSize + 3 * i + j;

              // Mark for removal
              if (!keepBuilding[marker])
                keepCell[cellIndex] = false;
            }
          }
        }
      }

      // Phase 2: Extract new mesh for all cells that should be kept

      // Renumber vertices and cells
      std::unordered_map<size_t, size_t> vertexMap;
      std::unordered_map<size_t, size_t> cellMap;
      size_t k = 0;
      size_t l = 0;
      for (size_t i = 0; i < mesh3D.Cells.size(); i++)
      {
        if (keepCell[i])
        {
          // Get cell
          const Simplex3D &T = mesh3D.Cells[i];

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
          cellMap[i] = l++;
        }
      }

      // Initialize new mesh data
      const size_t numVertices = vertexMap.size();
      const size_t numCells = cellMap.size();
      std::vector<Point3D> vertices(numVertices);
      std::vector<Simplex3D> cells(numCells);
      std::vector<int> markers(numCells);

      // Set new mesh data
      for (const auto v : vertexMap)
        vertices[v.second] = mesh3D.Vertices[v.first];
      for (const auto c : cellMap)
      {
        cells[c.second].v0 = vertexMap[mesh3D.Cells[c.first].v0];
        cells[c.second].v1 = vertexMap[mesh3D.Cells[c.first].v1];
        cells[c.second].v2 = vertexMap[mesh3D.Cells[c.first].v2];
        cells[c.second].v3 = vertexMap[mesh3D.Cells[c.first].v3];
        markers[c.second] = mesh3D.Markers[c.first];
      }

      // Assign new mesh data
      mesh3D.Vertices = vertices;
      mesh3D.Cells = cells;
      mesh3D.Markers = markers;
    }

    // Generate 3D surface meshes for visualization. The first surface is
    // the ground (height map) and the remaining surfaces are the extruded
    // building footprints. Note that meshes are non-conforming; the ground
    // and building meshes are non-matching and the building meshes may
    // contain hanging nodes.
    static std::vector<Surface3D> GenerateSurfaces3D(const CityModel &cityModel,
                                                     const GridField2D &dtm,
                                                     double xMin,
                                                     double yMin,
                                                     double xMax,
                                                     double yMax,
                                                     double resolution,
                                                     bool flatGround)
    {
      Info("MeshGenerator: Generating 3D surface meshes...");
      Timer("GenerateSurfaces3D");

      // Create empty list of surfaces
      std::vector<Surface3D> surfaces;

      // Generate empty subdomains for Triangle mesh generation
      std::vector<std::vector<Point2D>> subDomains;

      // Generate boundary for Triangle mesh generation
      std::vector<Point2D> boundary;
      boundary.push_back(Point2D(xMin, yMin));
      boundary.push_back(Point2D(xMax, yMin));
      boundary.push_back(Point2D(xMax, yMax));
      boundary.push_back(Point2D(xMin, yMax));

      // Generate 2D mesh of domain
      Info("MeshGenerator: Generating ground mesh");
      Mesh2D mesh2D;
      CallTriangle(mesh2D, boundary, subDomains, resolution);

      // Compute domain markers
      ComputeDomainMarkers(mesh2D, cityModel);

      // Create ground surface with zero height
      Surface3D surface3D;
      surface3D.Faces = mesh2D.Cells;
      surface3D.Vertices.resize(mesh2D.Vertices.size());
      for (size_t i = 0; i < mesh2D.Vertices.size(); i++)
      {
        const Point2D &p2D = mesh2D.Vertices[i];
        Vector3D p3D(p2D.x, p2D.y, 0.0);
        surface3D.Vertices[i] = p3D;
      }

      // Displace ground surface
      Info("MeshGenerator: Displacing ground mesh");
      if (flatGround)
      {
        // If ground is flat, just iterate over vertices and set height
        const double z = dtm.Min();
        for (size_t i = 0; i < mesh2D.Vertices.size(); i++)
          surface3D.Vertices[i].z = z;
      }
      else
      {
        // Fill all points with maximum height. This is used to
        // always choose the smallest height for each point since
        // each point may be visited multiple times.
        const double zMax = dtm.Max();
        for (size_t i = 0; i < mesh2D.Vertices.size(); i++)
          surface3D.Vertices[i].z = zMax;

        // If ground is not float, iterate over the triangles
        for (size_t i = 0; i < mesh2D.Cells.size(); i++)
        {
          // Get cell marker
          const int cellMarker = mesh2D.Markers[i];

          // Get triangle
          const Simplex2D& T = mesh2D.Cells[i];

          // Check cell marker
          if (cellMarker != -2) // not ground
          {
            // Compute minimum height of vertices
            double zMin = std::numeric_limits<double>::max();
            zMin = std::min(zMin, dtm(mesh2D.Vertices[T.v0]));
            zMin = std::min(zMin, dtm(mesh2D.Vertices[T.v1]));
            zMin = std::min(zMin, dtm(mesh2D.Vertices[T.v2]));

            // Set minimum height for all vertices
            setMin(surface3D.Vertices[T.v0].z, zMin);
            setMin(surface3D.Vertices[T.v1].z, zMin);
            setMin(surface3D.Vertices[T.v2].z, zMin);
          }
          else
          {
            // Sample height map at vertex position for all vertices
            setMin(surface3D.Vertices[T.v0].z, dtm(mesh2D.Vertices[T.v0]));
            setMin(surface3D.Vertices[T.v1].z, dtm(mesh2D.Vertices[T.v1]));
            setMin(surface3D.Vertices[T.v2].z, dtm(mesh2D.Vertices[T.v2]));
          }
        }
      }

      // Add ground surface to array of surfaces
      surfaces.push_back(surface3D);

      // Get ground height (minimum)
      const double groundHeight = dtm.Min();

      // Iterate over buildings to generate surfaces
      Info("MeshGenerator: Generating building meshes");
      for (auto const &building : cityModel.Buildings)
      {
        // Generate 2D mesh of building footprint
        Mesh2D _mesh2D;
        CallTriangle(_mesh2D, building.Footprint.Vertices, subDomains,
                     resolution);

        // Create empty 3D surface
        Surface3D _surface3D;

        // Note: The generated 2D mesh contains all the input boundary
        // points with the same numbers as in the footprint polygon, but
        // may also contain new points (Steiner points) added during
        // mesh generation. We add the top points (including any Steiner
        // points) first, then the points at the bottom (the footprint).

        // Set height of building
        const double buildingHeight = building.MaxHeight();

        // Set total number of points
        const size_t numMeshPoints = _mesh2D.Vertices.size();
        const size_t numBoundaryPoints = building.Footprint.Vertices.size();
        _surface3D.Vertices.resize(numMeshPoints + numBoundaryPoints);

        // Set total number of triangles
        const size_t numMeshTriangles = _mesh2D.Cells.size();
        const size_t numBoundaryTriangles = 2 * numBoundaryPoints;
        _surface3D.Faces.resize(numMeshTriangles + numBoundaryTriangles);

        // Add points at top
        for (size_t i = 0; i < numMeshPoints; i++)
        {
          const Point2D &p2D = _mesh2D.Vertices[i];
          const Vector3D p3D(p2D.x, p2D.y, buildingHeight);
          _surface3D.Vertices[i] = p3D;
        }

        // Add points at bottom
        for (size_t i = 0; i < numBoundaryPoints; i++)
        {
          const Point2D &p2D = _mesh2D.Vertices[i];
          const Vector3D p3D(p2D.x, p2D.y, groundHeight);
          _surface3D.Vertices[numMeshPoints + i] = p3D;
        }

        // Add triangles on top
        for (size_t i = 0; i < numMeshTriangles; i++)
          _surface3D.Faces[i] = _mesh2D.Cells[i];

        // Add triangles on boundary
        for (size_t i = 0; i < numBoundaryPoints; i++)
        {
          const size_t v0 = i;
          const size_t v1 = (i + 1) % numBoundaryPoints;
          const size_t v2 = v0 + numMeshPoints;
          const size_t v3 = v1 + numMeshPoints;
          Simplex2D t0(v0, v2, v1); // Outward-pointing normal
          Simplex2D t1(v1, v2, v3); // Outward-pointing normal
          _surface3D.Faces[numMeshTriangles + 2 * i] = t0;
          _surface3D.Faces[numMeshTriangles + 2 * i + 1] = t1;
        }

        // Add surface
        surfaces.push_back(_surface3D);
      }

      return surfaces;
    }

  private:
    // Call Triangle to compute 2D mesh
    static void
    CallTriangle(Mesh2D &mesh2D,
                 const std::vector<Point2D> &boundary,
                 const std::vector<std::vector<Point2D>> &subDomains,
                 double h)
    {
      Timer("CallTriangle");

      // Set area constraint to control mesh size
      const double maxArea = 0.5 * h * h;

      // Set input switches for Triangle
      char triswitches[64];
      sprintf(triswitches, "zQpq25a%.16f", maxArea);

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
      //PrintTriangleIO(out);
      //PrintTriangleIO(vorout);

      // Extract points
      mesh2D.Vertices.reserve(out.numberofpoints);
      for (int i = 0; i < out.numberofpoints; i++)
      {
        Point2D p(out.pointlist[2 * i], out.pointlist[2 * i + 1]);
        mesh2D.Vertices.push_back(p);
      }

      // Extract triangles
      mesh2D.Cells.reserve(out.numberoftriangles);
      for (int i = 0; i < out.numberoftriangles; i++)
      {
        // Note the importance of creating a sorted simplex here!
        Simplex2D t(out.trianglelist[3 * i], out.trianglelist[3 * i + 1],
                    out.trianglelist[3 * i + 2], true);
        mesh2D.Cells.push_back(t);
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

      io.pointlist                  = nullptr;
      io.pointmarkerlist            = nullptr;
      io.pointmarkerlist            = nullptr;
      io.numberofpoints             = 0;
      io.numberofpointattributes    = 0;
      io.trianglelist               = nullptr;
      io.triangleattributelist      = nullptr;
      io.trianglearealist           = nullptr;
      io.neighborlist               = nullptr;
      io.numberoftriangles          = 0;
      io.numberofcorners            = 0;
      io.numberoftriangleattributes = 0;
      io.segmentlist                = nullptr;
      io.segmentmarkerlist          = nullptr;
      io.numberofsegments           = 0;
      io.holelist                   = nullptr;
      io.numberofholes              = 0;
      io.regionlist                 = nullptr;
      io.numberofregions            = 0;
      io.edgelist                   = nullptr;
      io.edgemarkerlist             = nullptr;
      io.normlist                   = nullptr;
      io.numberofedges              = 0;

      return io;
    }

    // Print triangle I/O data
    static void PrintTriangleIO(const struct triangulateio& io)
    {
      Info("Triangle I/O data: ");
      Info("  pointlist = " + str(reinterpret_cast<std::uintptr_t>(io.pointlist)));
      Info("  pointmarkerlist = " + str(reinterpret_cast<std::uintptr_t>(io.pointmarkerlist)));
      if (io.pointmarkerlist)
      {
        std::stringstream stringBuilder{};
        stringBuilder << "   ";
        for (int i = 0; i < io.numberofpoints; i++)
          stringBuilder << " " << io.pointmarkerlist[i];
        stringBuilder << std::endl;
        Info(stringBuilder.str());
      }
      Info("  numberofpoints = " + str(io.numberofpoints));
      Info("  numberofpointattributes = " + str(io.numberofpointattributes));
      Info("  trianglelist = " + str(reinterpret_cast<std::uintptr_t>(io.trianglelist)));
      Info("  triangleattributelist = " + str(reinterpret_cast<std::uintptr_t>(io.triangleattributelist)));
      Info("  trianglearealist = " + str(reinterpret_cast<std::uintptr_t>(io.trianglearealist)));
      Info("  neighborlist = " + str(reinterpret_cast<std::uintptr_t>(io.neighborlist)));
      Info("  numberoftriangles = " + str(io.numberoftriangles));
      Info("  numberofcorners = " + str(io.numberofcorners));
      Info("  numberoftriangleattributes = " + str(io.numberoftriangleattributes));
      Info("  segmentlist = " + str(reinterpret_cast<std::uintptr_t>(io.segmentlist)));
      Info("  segmentmarkerlist = " + str(reinterpret_cast<std::uintptr_t>(io.segmentmarkerlist)));
      if (io.segmentmarkerlist)
      {
        std::stringstream stringBuilder{};
        stringBuilder << "   ";
        for (int i = 0; i < io.numberofsegments; i++)
          stringBuilder << " " << io.segmentmarkerlist[i];
        stringBuilder << std::endl;
        Info(stringBuilder.str());
      }
      Info("  numberofsegments = " + str(io.numberofsegments));
      Info("  holelist = " + str(reinterpret_cast<std::uintptr_t>(io.holelist)));
      Info("  numberofholes = " + str(io.numberofholes));
      Info("  regionlist = " + str(reinterpret_cast<std::uintptr_t>(io.regionlist)));
      Info("  numberofregions = " + str(io.numberofregions));
      Info("  edgelist = " + str(reinterpret_cast<std::uintptr_t>(io.edgelist)));
      Info("  edgemarkerlist = " + str(reinterpret_cast<std::uintptr_t>(io.edgemarkerlist)));
      Info("  normlist = " + str(reinterpret_cast<std::uintptr_t>(io.normlist)));
      Info("  numberofedges = " + str(io.numberofedges));
    }

    // Compute domain markers for subdomains
    static void ComputeDomainMarkers(Mesh2D & mesh, const CityModel & cityModel)
    {
      Info("MeshGenerator: Computing domain markers");
      Timer timer("ComputeDomainMarkers");

      // Build search tree for city model
      cityModel.BuildSearchTree();

      // Initialize domain markers and set all markers to -2 (ground)
      mesh.Markers.resize(mesh.Cells.size());
      std::fill(mesh.Markers.begin(), mesh.Markers.end(), -2);

      // Initialize markers for vertices belonging to a building
      std::vector<bool> isBuildingVertex(mesh.Vertices.size());
      std::fill(isBuildingVertex.begin(), isBuildingVertex.end(), false);

      // Iterate over cells to mark buildings
      for (size_t i = 0; i < mesh.Cells.size(); i++)
      {
        // Find building containg midpoint of cell (if any)
        const Point2D c = mesh.MidPoint(i);
        const int marker = cityModel.FindBuilding(Vector2D(c));

        // Get triangle
        const Simplex2D &T = mesh.Cells[i];

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
        // (not only midpoint). Necessary for when generating
        // visualization meshes that are not boundary-fitted.
        if (cityModel.FindBuilding(Vector2D(mesh.Vertices[T.v0])) >= 0)
          isBuildingVertex[T.v0] = true;
        if (cityModel.FindBuilding(Vector2D(mesh.Vertices[T.v1])) >= 0)
          isBuildingVertex[T.v1] = true;
        if (cityModel.FindBuilding(Vector2D(mesh.Vertices[T.v2])) >= 0)
          isBuildingVertex[T.v2] = true;
      }

      // Iterate over cells to mark building halos
      for (size_t i = 0; i < mesh.Cells.size(); i++)
      {
        // Check if any of the cell vertices belongs to a building
        const Simplex2D &T = mesh.Cells[i];
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
    static void setMin(double& x, double y)
    {
      if (y < x)
        x = y;
    }
  };

} // namespace DTCC

#endif
