// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

// Developer note: The implementation of this class does not follow
// DTCC Coding Standard since the code makes heavy use of FEniCS/DOLFIN
// which uses a quite different coding standard and the mix looks horrible.

#ifndef DTCC_LAPLACIAN_SMOOTHER_H
#define DTCC_LAPLACIAN_SMOOTHER_H

#include <cmath>
#include <dolfin.h>
#include <iostream>

#include "Laplacian.h"
#include "LinearSpace2D.h"
#include "Logging.h"
#include "Mesh.h"
#include "Timer.h"

namespace DTCC
{

  class LaplacianSmoother
  {
  public:

    // Smooth mesh using Laplacian smoothing
    static void SmoothMesh(dolfin::Mesh &mesh,
                           const GridField2D &heightMap,
                           const CityModel &cityModel,
                           const std::vector<int> &domainMarkers,
                           double h)
    {
      Info("Smoothing mesh (Laplacian smoothing)...");
      Timer("SmoothMesh");

      // Get number of vertices
      const size_t num_vertices = mesh.num_vertices();

      // Create function space and bilinear form
      auto m = std::make_shared<dolfin::Mesh>(mesh);
      auto V = std::make_shared<Laplacian::FunctionSpace>(m);
      auto a = std::make_shared<Laplacian::BilinearForm>(V, V);

      // Assemble matrix
      auto A = std::make_shared<dolfin::Matrix>();
      dolfin::assemble(*A, *a);

      // Initialize vectors
      auto x = std::make_shared<dolfin::Vector>();
      auto b = std::make_shared<dolfin::Vector>();
      A->init_vector(*x, 0);
      A->init_vector(*b, 0);

      // Create boundary markers from domain markers
      auto subDomains = std::make_shared<dolfin::MeshFunction<size_t>>(m, 2);
      ComputeBoundaryMarkers(*subDomains, domainMarkers);

      // Create expressions for boundary values (heights)
      auto h0 = std::make_shared<GroundExpression>(heightMap);
      auto h1 = std::make_shared<HaloExpression>(heightMap, mesh);
      auto h2 = std::make_shared<BuildingsExpression>(cityModel, domainMarkers);

      // Create boundary conditions
      auto bc0 = std::make_shared<dolfin::DirichletBC>(V, h0, subDomains, 2);
      auto bc1 = std::make_shared<dolfin::DirichletBC>(V, h1, subDomains, 0);
      auto bc2 = std::make_shared<dolfin::DirichletBC>(V, h2, subDomains, 0);

      // Apply boundary conditions
      bc2->apply(*A, *b);
      bc1->apply(*A, *b);
      bc0->apply(*A, *b);

      // Create linear solver
      dolfin::KrylovSolver solver(mesh.mpi_comm(), "bicgstab", "amg");
      solver.parameters["nonzero_initial_guess"] = true;
      solver.set_operator(A);

      // Solve linear system
      *x = *b;
      solver.solve(*x, *b);

      // Get displacement
      const std::vector<dolfin::la_index> v2d = vertex_to_dof_map(*V);
      std::vector<double> dz(num_vertices);
      x->get_local(dz.data(), num_vertices, v2d.data());

      // Update mesh coordinates
      double coordinates[3];
      for (std::size_t i = 0; i < num_vertices; i++)
      {
        coordinates[0] = mesh.geometry().x(i, 0);
        coordinates[1] = mesh.geometry().x(i, 1);
        coordinates[2] = mesh.geometry().x(i, 2) + dz[i];
        mesh.geometry().set(i, coordinates);
      }
    }

    // Smooth mesh using elastic smoothing
    static void SmoothMeshElastic(dolfin::Mesh &mesh,
                                  const GridField2D &heightMap,
                                  const CityModel &cityModel,
                                  const std::vector<int> &domainMarkers,
                                  double h)
    {
      std::cout << "Elastic smoothing not (yet) implemented." << std::endl;
    }

    // Generate height map function (used only for testing/visualization)
    static std::shared_ptr<dolfin::Function>
    GenerateHeightMapFunction(const dolfin::Mesh &mesh,
                              const GridField2D &heightMap)
    {
      // Create function space
      auto m = std::make_shared<dolfin::Mesh>(mesh);
      auto V = std::make_shared<LinearSpace2D::FunctionSpace>(m);

      // Create boundary condition
      auto bcz = std::make_shared<dolfin::DirichletBC>(
                                                       V, std::make_shared<GroundExpression>(heightMap),
                                                       std::make_shared<EntireDomain>());

      // Create function and apply boundary condition
      auto z = std::make_shared<dolfin::Function>(V);
      bcz->apply(*z->vector());

      return z;
    }

  private:
    // Tolerance for geometric tests
    static constexpr double tol = 1e-3;

    // Boundary definition for entire domain
    class EntireDomain : public dolfin::SubDomain
    {
      bool inside(const dolfin::Array<double> &x, bool on_boundary) const
      {
        return true;
      }
    };

    // A note on subtracting z-coordinate in Expressions below: Since
    // we solve for the z-displacement, we need to subtract the z-coordinate
    // but since the expressions are also used by GenerateHeightMapFunction()
    // to generate a 2D height map with absolute height values, the
    // subtraction should only be done in the 3D case.

    // Boundary value for buildings
    class BuildingsExpression : public dolfin::Expression
    {
    public:
      // Reference to city model
      const CityModel &cityModel;

      // Reference to domain markers
      const std::vector<int> &domainMarkers;

      // Building heights (absolute z-coordinates of roofs)
      std::vector<double> buildingHeights;

      // Constructor
      BuildingsExpression(const CityModel &cityModel,
                          const std::vector<int> &domainMarkers)
        : cityModel(cityModel), domainMarkers(domainMarkers)
      {
      }

      // Evaluation of z-displacement
      void eval(dolfin::Array<double> &values,
                const dolfin::Array<double> &x,
                const ufc::cell &ufc_cell) const
      {
        // Get building height
        const size_t i = domainMarkers[ufc_cell.index];
        const double z = cityModel.Buildings[i].Height;

        // Set height of building
        values[0] = z;

        // See note above on subtracting z-coordinate
        if (x.size() == 3)
          values[0] -= x[2];
      }
    };

    // Boundary value for building halos (2D)
    class HaloExpression : public dolfin::Expression
    {
    public:
      // Reference to height map
      const GridField2D &heightMap;

      // Reference to mesh
      const dolfin::Mesh &mesh;

      // Constructor
      HaloExpression(const GridField2D &heightMap, const dolfin::Mesh &mesh)
        : Expression(), heightMap(heightMap), mesh(mesh)
      {
      }

      // Evaluation of z-displacement
      void eval(dolfin::Array<double> &values,
                const dolfin::Array<double> &x,
                const ufc::cell &ufc_cell) const
      {
        // Get DOLFIN cell
        dolfin::Cell cell(mesh, ufc_cell.index);

        // Check each cell vertex and take the minimum
        double zmin = std::numeric_limits<double>::max();
        for (dolfin::VertexIterator vertex(cell); !vertex.end(); ++vertex)
        {
          const dolfin::Point p = vertex->point();
          const double z = heightMap(Vector2D(p.x(), p.y()));
          zmin = std::min(zmin, z);
        }

        // Set value
        values[0] = zmin;

        // See note above on subtracting z-coordinate
        if (x.size() == 3)
          values[0] -= x[2];
      }
    };

    // Boundary value for ground
    class GroundExpression : public dolfin::Expression
    {
    public:
      // Reference to height map
      const GridField2D &heightMap;

      // Constructor
      GroundExpression(const GridField2D &heightMap)
        : Expression(), heightMap(heightMap)
      {
      }

      // Evaluation of z-displacement
      void eval(dolfin::Array<double> &values,
                const dolfin::Array<double> &x) const
      {
        // Evaluate height map
        values[0] = heightMap(Vector2D(x[0], x[1]));

        // See note above on subtracting z-coordinate
        if (x.size() == 3)
          values[0] -= x[2];
      }
    };

    // Compute boundary markers from domain markers
    static void ComputeBoundaryMarkers(dolfin::MeshFunction<size_t> &subDomains,
                                       std::vector<int> domainMarkers)
    {
      // The domain markers indicate a nonnegative building number for cells
      // that touch the roofs of buildings, -1 for cells that touch the
      // ground close to buildings, -2 for other cells that touch the ground
      // and -3 for remaining cells. We now need to convert these *cell*
      // markers to *facet* markers. The facet markers are set as follows:
      //
      // 0: roofs of buildings
      // 1: ground close to buildings (converted from -1 cell markers)
      // 2: ground away from buildings (converted from -2 cell markers)
      // 3: everything else

      // Iterate over the facets of the mesh
      for (dolfin::FacetIterator f(*subDomains.mesh()); !f.end(); ++f)
      {
        // Set default facet marker
        size_t facetMarker = 3;

        // Check if we are on the boundary
        if (f->exterior())
        {
          // Get index of neighboring cell (should only be one)
          const size_t cellIndex = f->entities(3)[0];

          // Check z-component of cell normal
          const double nz = f->normal(2);
          const bool downwardFacet = nz <= -1.0 + tol;

          // Get cell marker
          const int cellMarker = domainMarkers[cellIndex];

          // Set facet marker (default to 3 in the else case)
          if (downwardFacet && cellMarker >= 0)
            facetMarker = 0;
          else if (downwardFacet && cellMarker == -1)
            facetMarker = 1;
          else if (downwardFacet && cellMarker == -2)
            facetMarker = 2;
        }

        // Set marker value
        subDomains.set_value(f->index(), facetMarker);
      }
    }
  };

} // namespace DTCC

#endif
