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
                           double H, double h)
    {
        SmoothMeshLaplacian(mesh, heightMap, cityModel, domainMarkers, H, h);
    }

    // Smooth mesh using Laplacian smoothing
    static void SmoothMeshLaplacian(dolfin::Mesh& mesh,
                                    const HeightMap& heightMap,
                                    const CityModel& cityModel,
                                    const std::vector<int>& domainMarkers,
                                    double H, double h)
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

        std::make_shared<BuildingsDomain>(H, h);

        // Create boundary conditions
        auto bcg = std::make_shared<dolfin::DirichletBC>
                   (V,
                    std::make_shared<GroundExpression>(heightMap),
                    std::make_shared<GroundDomain>());
        auto bcb = std::make_shared<dolfin::DirichletBC>
                   (V,
                    std::make_shared<BuildingsExpression>(heightMap,
                            cityModel,
                            domainMarkers),
                    std::make_shared<BuildingsDomain>(H, h));

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

        // Get new z-coordinates
        const std::vector<dolfin::la_index> v2d = vertex_to_dof_map(*V);
        std::vector<double> z(num_vertices);
        x->get_local(z.data(), num_vertices, v2d.data());

        // Update mesh coordinates
        double coordinates[3];
        for (std::size_t i = 0; i < num_vertices; i++)
        {
            coordinates[0] = mesh.geometry().x(i, 0);
            coordinates[1] = mesh.geometry().x(i, 1);
            coordinates[2] = z[i];
            mesh.geometry().set(i, coordinates);
        }
    }

    // Smooth mesh using elastic smoothing
    static void SmoothMeshElastic(dolfin::Mesh& mesh,
                                  const HeightMap& heightMap,
                                  const CityModel& cityModel,
                                  const std::vector<int>& domainMarkers,
                                  double H, double h)
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
                    std::make_shared<GroundExpression>(heightMap),
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

    // Boundary definition for ground (height map)
    class GroundDomain : public dolfin::SubDomain
    {
        bool inside(const dolfin::Array<double>& x, bool on_boundary) const
        {
            return on_boundary && x[2] < tol;
        }
    };

    // Boundary definition for buildings
    class BuildingsDomain : public dolfin::SubDomain
    {
    public:

        // Domain height and mesh size
        double H, h;

        // Constructor
        BuildingsDomain(double H, double h) : H(H), h(h) {}

        // We use a "clever" trick to specify building roofs geometrically
        bool inside(const dolfin::Array<double>& x, bool on_boundary) const
        {
            const bool onGrid = std::fmod(x[2], h) < tol;
            const bool onBottom = x[2] < tol;
            const bool onTop = x[2] > H - tol;
            return on_boundary && onGrid && !onBottom && !onTop;
        }

    };

    // Boundary value for ground (height map)
    class GroundExpression : public dolfin::Expression
    {
    public:

        // Reference to height map
        const HeightMap& heightMap;

        // Constructor
        GroundExpression(const HeightMap& heightMap)
            : heightMap(heightMap), Expression() {}

        // Evaluation
        void eval(dolfin::Array<double>& values,
                  const dolfin::Array<double>& x) const
        {
            values[0] = heightMap(x[0], x[1]);
        }

    };

    // Boundary value for buildings
    class BuildingsExpression : public dolfin::Expression
    {
    public:

        // Reference to domain markers
        const std::vector<int>& domainMarkers;

        // Building heights (absolute z-coordinates of roofs)
        std::vector<double> buildingHeights;

        // Constructor
        BuildingsExpression(const HeightMap& heightMap,
                            const CityModel& cityModel,
                            const std::vector<int>& domainMarkers)
            : domainMarkers(domainMarkers),
              buildingHeights(cityModel.Buildings.size()),
              Expression()
        {
            // Compute height of each building
            for (size_t i = 0; i < cityModel.Buildings.size(); i++)
            {
                // Sample height map at center
                const Point2D c = cityModel.Buildings[i].Center();
                const double z0 = heightMap(c.x, c.y);

                // Add relative building height
                buildingHeights[i] = z0 + 30.0; //cityModel.Buildings.Height;
            }
        }

        // Evaluation
        void eval(dolfin::Array<double>& values,
                  const dolfin::Array<double>& x,
                  const ufc::cell& cell) const
        {
            // Get number of building
            const int i = domainMarkers[cell.index];

            // FIXME: Why do we get -2?

            const double z = (i >= 0 ? buildingHeights[i] : 0.0);

            // Set height of building
            values[0] = z;
        }

    };

};

}

#endif
