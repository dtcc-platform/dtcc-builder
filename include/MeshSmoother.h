// PDE based mesh smoothing.
// Copyright (C) 2018 Anders Logg.

// Developer note: The implementation of this class does not follow
// UE4 Coding Standard since the code makes heavy use of FEniCS/DOLFIN
// which uses a quite different coding standard and the mix looks horrible.

#ifndef VC_MESH_SMOOTHER_H
#define VC_MESH_SMOOTHER_H

#include <iostream>
#include <cmath>
#include <dolfin.h>

#include "Mesh.h"
#include "HeightMap.h"
#include "LaplacianSmoother.h"
#include "LinearSpace2D.h"

namespace VirtualCity
{

class MeshSmoother
{
public:

    // Smooth mesh using default method (Laplacian smoothing)
    static void SmoothMesh(dolfin::Mesh& mesh,
                           const HeightMap& heightMap,
                           const CityModel& cityModel,
                           const std::vector<int>& domainMarkers,
                           double h)
    {
        SmoothMeshLaplacian(mesh, heightMap, cityModel, domainMarkers, h);
    }

    // Smooth mesh using Laplacian smoothing
    static void SmoothMeshLaplacian(dolfin::Mesh& mesh,
                                    const HeightMap& heightMap,
                                    const CityModel& cityModel,
                                    const std::vector<int>& domainMarkers,
                                    double h)
    {
        std::cout << "Smoothing mesh (Laplacian smoothing)..." << std::endl;

        // Get mesh sizes
        const size_t num_cells = mesh.num_cells();
        const size_t num_vertices = mesh.num_vertices();

        // Create function space and bilinear form
        auto m = std::make_shared<dolfin::Mesh>(mesh);
        auto V = std::make_shared<LaplacianSmoother::FunctionSpace>(m);
        auto a = std::make_shared<LaplacianSmoother::BilinearForm>(V, V);

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

        // Create expressions for ground and building heights
        auto hg = std::make_shared<GroundExpression3D>(heightMap);
        auto hb = std::make_shared<BuildingsExpression>
                  (cityModel, domainMarkers);

        // Create boundary conditions
        auto bcg = std::make_shared<dolfin::DirichletBC>
                   (V, hg, subDomains, 1);
        auto bcb = std::make_shared<dolfin::DirichletBC>
                   (V, hb, subDomains, 0);

        // Apply boundary conditions
        bcg->apply(*A, *b);
        bcb->apply(*A, *b);

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
    static void SmoothMeshElastic(dolfin::Mesh& mesh,
                                  const HeightMap& heightMap,
                                  const CityModel& cityModel,
                                  const std::vector<int>& domainMarkers,
                                  double h)
    {
        std::cout << "Elastic smoothing not (yet) implemented." << std::endl;
    }

    // Generate height map function (used only for testing/visualization)
    static std::shared_ptr<dolfin::Function>
    GenerateHeightMapFunction(const dolfin::Mesh& mesh,
                              const HeightMap& heightMap)
    {
        // Create function space
        auto m = std::make_shared<dolfin::Mesh>(mesh);
        auto V = std::make_shared<LinearSpace2D::FunctionSpace>(m);

        // Create boundary condition
        auto bcz = std::make_shared<dolfin::DirichletBC>
                   (V,
                    std::make_shared<GroundExpression2D>(heightMap),
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
        bool inside(const dolfin::Array<double>& x, bool on_boundary) const
        {
            return true;
        }
    };

    // Boundary value for ground (2D)
    class GroundExpression2D : public dolfin::Expression
    {
    public:

        // Reference to height map
        const HeightMap& heightMap;

        // Constructor
        GroundExpression2D(const HeightMap& heightMap)
            : heightMap(heightMap), Expression() {}

        // Evaluation of z-displacement
        void eval(dolfin::Array<double>& values,
                  const dolfin::Array<double>& x) const
        {
            values[0] = heightMap(x[0], x[1]);
        }

    };

    // Boundary value for ground (3D)
    class GroundExpression3D : public dolfin::Expression
    {
    public:

        // Reference to height map
        const HeightMap& heightMap;

        // Constructor
        GroundExpression3D(const HeightMap& heightMap)
            : heightMap(heightMap), Expression() {}

        // Evaluation of z-displacement
        void eval(dolfin::Array<double>& values,
                  const dolfin::Array<double>& x) const
        {
            values[0] = heightMap(x[0], x[1]) - x[2];
        }

    };

    // Boundary value for buildings
    class BuildingsExpression : public dolfin::Expression
    {
    public:

        // Reference to city model
        const CityModel& cityModel;

        // Reference to domain markers
        const std::vector<int>& domainMarkers;

        // Building heights (absolute z-coordinates of roofs)
        std::vector<double> buildingHeights;

        // Constructor
        BuildingsExpression(const CityModel& cityModel,
                            const std::vector<int>& domainMarkers)
            : cityModel(cityModel), domainMarkers(domainMarkers) {}

        // Evaluation of z-displacement
        void eval(dolfin::Array<double>& values,
                  const dolfin::Array<double>& x,
                  const ufc::cell& cell) const
        {
            // Get building height
            const size_t i = domainMarkers[cell.index];
            const double z = cityModel.Buildings[i].Height;

            // Set height of building
            values[0] = z - x[2];
        }

    };

    // Compute boundary markers from domain markers
    static void
    ComputeBoundaryMarkers(dolfin::MeshFunction<size_t>& subDomains,
                           std::vector<int> domainMarkers)
    {
        // The domain markers indicate a nonnegative building number for cells
        // that touch the roofs of buildings, -1 for cells that touch the
        // ground, and -2 for other cells. We now need to convert these *cell*
        // markers to *facet* markers.
        //
        // Facets that touch the roofs of buildings are marked as 0, facets
        // touching the ground are marked as 1 and the rest are marked as 2.

        // Iterate over the facets of the mesh
        for (dolfin::FacetIterator f(*subDomains.mesh()); !f.end(); ++f)
        {
            // Set default facet marker
            size_t facetMarker = 2;

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

                // Set facet marker (default to 2 in the else case)
                if (downwardFacet && cellMarker >= 0)
                    facetMarker = 0;
                else if (downwardFacet && cellMarker == -1)
                    facetMarker = 1;
            }

            // Set marker value
            subDomains.set_value(f->index(), facetMarker);
        }
    }

};

}

#endif
