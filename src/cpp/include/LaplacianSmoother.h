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

#include "FEniCS.h"
#include "Laplacian.h"
#include "LinearSpace2D.h"
#include "Logging.h"
#include "Mesh.h"
#include "Timer.h"

namespace DTCC_BUILDER
{

class LaplacianSmoother
{
public:
  // Smooth mesh using Laplacian smoothing
  static void SmoothMesh3D(Mesh3D &mesh3D,
                           const CityModel &city_model,
                           const GridField2D &dem,
                           double topHeight,
                           bool fixBuildings,
                           bool writeMatrix)
  {
    info("LaplacianSmoother: Smoothing mesh (Laplacian smoothing)...");
    Timer timer("SmoothMesh3D");

    // Convert to FEniCS mesh
    dolfin::Mesh _mesh3D;
    FEniCS::ConvertMesh(mesh3D, _mesh3D);

    // Get number of vertices
    const size_t num_vertices = _mesh3D.num_vertices();

    // Create function space and bilinear form
    auto m = std::make_shared<dolfin::Mesh>(_mesh3D);
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
    ComputeBoundaryMarkers(*subDomains, mesh3D.Markers);

    // Create expressions for boundary values (heights)
    const auto h0 =
        std::make_shared<BuildingsExpression>(city_model, mesh3D.Markers);
    const auto h1 = std::make_shared<HaloExpressionDEM>(dem, _mesh3D);
    const auto h2 = std::make_shared<GroundExpression>(dem);
    const auto h3 = std::make_shared<TopExpression>(topHeight);

    // Create boundary conditions
    const auto bc0 =
        std::make_shared<dolfin::DirichletBC>(V, h0, subDomains, 0);
    const auto bc1 =
        std::make_shared<dolfin::DirichletBC>(V, h1, subDomains, 1);
    const auto bc2 =
        std::make_shared<dolfin::DirichletBC>(V, h2, subDomains, 2);
    const auto bc3 =
        std::make_shared<dolfin::DirichletBC>(V, h3, subDomains, 3);

    // Note order of boundary conditions!

    // Apply boundary condition for buildings
    if (fixBuildings)
      bc0->apply(*A, *b);

    // Apply boundary conditions for top, ground and halos
    bc3->apply(*A, *b);
    bc2->apply(*A, *b);
    bc1->apply(*A, *b);

    // Write matrix to file (used for testing/development)
    // Loadin data in MATLAB:
    //
    // load matrix.dat
    // load vector.dat
    // A = spconvert(matrix);
    // b = vector;
    if (writeMatrix)
    {
      // Write matrix
      {
        const size_t M = A->size(0);
        const size_t N = A->size(1);
        info("LaplacianSmoother: Writing matrix of size " + str(M) + " x " +
             str(N) + " to file matrix.dat");
        std::stringstream ss;
        for (size_t i = 0; i < M; i++)
        {
          std::vector<std::size_t> columns{};
          std::vector<double> values{};
          A->getrow(i, columns, values);
          for (size_t j = 0; j < columns.size(); j++)
            ss << (i + 1) << " " << (columns[j] + 1) << " " << values[j]
               << std::endl;
        }
        std::ofstream fs;
        fs.open("matrix.dat");
        fs << ss.rdbuf();
        fs.close();
      }

      // Write vector
      {
        const size_t M = A->size(0);
        info("LaplacianSmoother: Writing vector of size " + str(M) +
             " to file vector.dat");
        std::stringstream ss;
        std::vector<double> values{};
        b->get_local(values);
        for (size_t i = 0; i < M; i++)
          ss << values[i] << std::endl;
        std::ofstream fs;
        fs.open("vector.dat");
        fs << ss.rdbuf();
        fs.close();
      }
    }

    // Create linear solver
    dolfin::KrylovSolver solver(_mesh3D.mpi_comm(), "bicgstab", "amg");
    solver.parameters["nonzero_initial_guess"] = true;
    solver.parameters["monitor_convergence"] = true;
    solver.set_operator(A);

    // Solve linear system
    *x = *b;
    solver.solve(*x, *b);

    // Get displacement
    const std::vector<dolfin::la_index> v2d = vertex_to_dof_map(*V);
    std::vector<double> dz(num_vertices);
    x->get_local(dz.data(), num_vertices, v2d.data());

    // Update mesh coordinates
    for (std::size_t i = 0; i < mesh3D.Vertices.size(); i++)
      mesh3D.Vertices[i].z += dz[i];
  }

  // Smooth mesh using elastic smoothing
  static void SmoothMeshElastic(dolfin::Mesh &mesh,
                                const GridField2D &dem,
                                const CityModel &city_model,
                                const std::vector<int> &domainMarkers,
                                double h)
  {
    warning("Elastic smoothing not (yet) implemented.");
  }

  // Generate elevation function (used only for testing/visualization)
  static std::shared_ptr<dolfin::Function>
  GenerateElevationFunction(const dolfin::Mesh &mesh,
                            const GridField2D &elevation)
  {
    // Create function space
    auto m = std::make_shared<dolfin::Mesh>(mesh);
    auto V = std::make_shared<LinearSpace2D::FunctionSpace>(m);

    // Create boundary condition
    auto bcz = std::make_shared<dolfin::DirichletBC>(
        V, std::make_shared<GroundExpression>(elevation),
        std::make_shared<EntireDomain>());

    // Create function and apply boundary condition
    auto z = std::make_shared<dolfin::Function>(V);
    bcz->apply(*z->vector());

    return z;
  }

private:
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
  // but since the expressions are also used by GenerateElevationFunction()
  // to generate a 2D height map with absolute height values, the
  // subtraction should only be done in the 3D case.

  // Boundary value for buildings
  class BuildingsExpression : public dolfin::Expression
  {
  public:
    // Reference to city model
    const CityModel &city_model;

    // Reference to domain markers
    const std::vector<int> &domainMarkers;

    // Building heights (absolute z-coordinates of roofs)
    std::vector<double> buildingHeights;

    // Constructor
    BuildingsExpression(const CityModel &city_model,
                        const std::vector<int> &domainMarkers)
        : city_model(city_model), domainMarkers(domainMarkers)
    {
    }

    // Evaluation of z-displacement
    void eval(dolfin::Array<double> &values,
              const dolfin::Array<double> &x,
              const ufc::cell &ufc_cell) const
    {
      // Get building height
      const size_t i = domainMarkers[ufc_cell.index];
      const double z = city_model.Buildings[i].MaxHeight();

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
    // Reference to city model
    const CityModel &city_model;

    // Reference to domain markers
    const std::vector<int> &domainMarkers;

    // Building heights (absolute z-coordinates of roofs)
    std::vector<double> buildingHeights;

    // Constructor
    HaloExpression(const CityModel &city_model,
                   const std::vector<int> &domainMarkers)
        : city_model(city_model), domainMarkers(domainMarkers)
    {
    }

    // Evaluation of z-displacement
    void eval(dolfin::Array<double> &values,
              const dolfin::Array<double> &x,
              const ufc::cell &ufc_cell) const
    {
      // Get building height
      const size_t i =
          domainMarkers[ufc_cell.index]; // FIXME: How to get building index?
      const double z = city_model.Buildings[i].MinHeight();

      // Set height of building
      values[0] = z;

      // See note above on subtracting z-coordinate
      if (x.size() == 3)
        values[0] -= x[2];
    }
  };

  // Boundary value for building halos based on DEM (2D)
  class HaloExpressionDEM : public dolfin::Expression
  {
  public:
    // Reference to elevation
    const GridField2D &elevation;

    // Reference to mesh
    const dolfin::Mesh &mesh;

    // Constructor
    HaloExpressionDEM(const GridField2D &elevation, const dolfin::Mesh &mesh)
        : Expression(), elevation(elevation), mesh(mesh)
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
        const double z = elevation(Vector2D(p.x(), p.y()));
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
    // Reference to elevation
    const GridField2D &elevation;

    // Constructor
    GroundExpression(const GridField2D &elevation)
        : Expression(), elevation(elevation)
    {
    }

    // Evaluation of z-displacement
    void eval(dolfin::Array<double> &values,
              const dolfin::Array<double> &x) const
    {
      // Evaluate height map
      values[0] = elevation(Vector2D(x[0], x[1]));

      // See note above on subtracting z-coordinate
      if (x.size() == 3)
        values[0] -= x[2];
    }
  };

  // Boundary value for top (no displacement)
  class TopExpression : public dolfin::Expression
  {
  public:
    // Absolute level of top
    double topHeight{};

    // Constructor
    TopExpression(double topHeight) : Expression(), topHeight(topHeight) {}

    // Evaluation of z-displacement
    void eval(dolfin::Array<double> &values,
              const dolfin::Array<double> &x) const
    {
      // Set absolute level of top
      values[0] = topHeight;

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
    // ground close to buildings, -2 for other cells that touch the ground,
    // -3 for top layer and -4 for remaining cells. We now need to convert
    // these *cell* markers to *facet* markers. The facet markers are set
    // as follows:
    //
    // 0: roofs of buildings (converted from non-negative cell markers)
    // 1: building halos (converted from -1 cell markers)
    // 2: ground (converted from -2 cell markers)
    // 3: top of domain (converted from -3 cell markers)
    // 4: everything else

    // Counters for boundary facets
    size_t k0 = 0;
    size_t k1 = 0;
    size_t k2 = 0;
    size_t k3 = 0;

    // Tolerance for checking whether we are on the top or bottom boundary.
    // Since we only displace the z-coordinates, the faces on the side should
    // always have a z-coordinate that is exactly zero.
    const double tol = 1e-6;

    // Iterate over the facets of the mesh
    for (dolfin::FacetIterator f(*subDomains.mesh()); !f.end(); ++f)
    {
      // Set default facet marker
      size_t facetMarker = 4;

      // Check if we are on the top or bottom boundary
      if (f->exterior() && std::abs(f->normal(2)) > tol)
      {
        // Get index of neighboring cell (should only be one)
        const size_t cellIndex = f->entities(3)[0];

        // Get cell marker
        const int cellMarker = domainMarkers[cellIndex];

        // Set facet marker
        if (cellMarker >= 0)
        {
          facetMarker = 0; // building
          k0++;
        }
        else if (cellMarker == -1)
        {
          facetMarker = 1; // halo
          k1++;
        }
        else if (cellMarker == -2)
        {
          facetMarker = 2; // ground
          k2++;
        }
        else if (cellMarker == -3)
        {
          facetMarker = 3; // top
          k3++;
        }
      }

      // Set marker value
      subDomains.set_value(f->index(), facetMarker);
    }

    const size_t k4 = subDomains.size() - (k0 + k1 + k2 + k3);
    info("LaplacianSmoother: Found " + str(k0) + " building boundary facets");
    info("LaplacianSmoother: Found " + str(k1) + " ground boundary facets");
    info("LaplacianSmoother: Found " + str(k2) + " halo boundary facets");
    info("LaplacianSmoother: Found " + str(k3) + " top boundary facets");
    info("LaplacianSmoother: Found " + str(k4) + " Neumann boundary facets");
  }
};

} // namespace DTCC_BUILDER

#endif
