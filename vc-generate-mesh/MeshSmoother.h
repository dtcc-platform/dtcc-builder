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
                           const std::vector<int>& domainMarkers)
    {
        SmoothMeshLaplacian(mesh, heightMap, domainMarkers);
    }

    // Smooth mesh using Laplacian smoothing
    static void SmoothMeshLaplacian(dolfin::Mesh& mesh,
                                    const HeightMap& heightMap,
                                    const std::vector<int>& domainMarkers)
    {
        std::cout << "Smoothing mesh (Laplacian smoothing)..." << std::endl;

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

        // Create boundary condition
        auto bcz = std::make_shared<dolfin::DirichletBC>
                   (V,
                    std::make_shared<HeightMapExpression>(heightMap),
                    std::make_shared<Bottom>());

        // Apply boundary conditions
        bcz->apply(*A, *b);

        // Create linear solver
        dolfin::KrylovSolver solver(mesh.mpi_comm(), "bicgstab", "amg");
        solver.parameters["nonzero_initial_guess"] = true;
        solver.set_operator(A);

        // Solve linear system
        *x = *b;
        solver.solve(*x, *b);

        // Get coordinate displacements
        const std::vector<dolfin::la_index> v2d = vertex_to_dof_map(*V);
        const size_t num_vertices = mesh.num_vertices();
        std::vector<double> displacements(num_vertices);
        x->get_local(displacements.data(), num_vertices, v2d.data());

        // Update mesh coordinates
        double coordinates[3];
        for (std::size_t i = 0; i < num_vertices; i++)
        {
            coordinates[0] = mesh.geometry().x(i, 0);
            coordinates[1] = mesh.geometry().x(i, 1);
            coordinates[2] = mesh.geometry().x(i, 2) + displacements[i];
            mesh.geometry().set(i, coordinates);
        }
    }

    // Smooth mesh using elastic smoothing
    static void SmoothMeshElastic(dolfin::Mesh& mesh,
                                  const HeightMap& heightMap,
                                  const std::vector<int>& domainMarkers)
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
                    std::make_shared<HeightMapExpression>(heightMap),
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
        bool inside(const dolfin::Array<double>& x, bool on_boundary) const
        {
            return true;
        }
    };

    // Boundary definition for bottom
    class Bottom : public dolfin::SubDomain
    {
        bool inside(const dolfin::Array<double>& x, bool on_boundary) const
        {
            return std::abs(x[2] - 0.0) < DOLFIN_EPS;
        }
    };

    // Boundary value for height map
    class HeightMapExpression : public dolfin::Expression
    {
    public:

        // Reference to actual height map
        const HeightMap& heightMap;

        // Create height map expression
        HeightMapExpression(const HeightMap& heightMap)
            : heightMap(heightMap), Expression()
        {
            // Do nothing
        }

        // Evaluation of height map
        void eval(dolfin::Array<double>& values,
                  const dolfin::Array<double>& x) const
        {
            values[0] = heightMap(x[0], x[1]);
        }

    };

};

}

#endif
