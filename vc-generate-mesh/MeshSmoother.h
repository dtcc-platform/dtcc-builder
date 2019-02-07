// PDE based mesh smoothing
// Copyright (C) 2018 Anders Logg.

#ifndef MESH_SMOOTHER_H
#define MESH_SMOOTHER_H

#include <iostream>
#include <dolfin.h>

#include "LaplacianSmoother.h"

namespace VirtualCity
{

class MeshSmoother
{
public:

    // Smooth mesh using default method (Laplacian smoothing)
    static void SmoothMesh()
    {
        SmoothMeshLaplacian();
    }

    // Smooth mesh using Laplacian smoothing
    static void SmoothMeshLaplacian()
    {
        std::cout << "Smoothing mesh (Laplacian smoothing)..." << std::endl;

        // FIXME: Test data
        //-----------------------------------------------------
        dolfin::Mesh mesh("Mesh2D.xml");




    }

    // Smooth mesh using elastic smoothing
    static void SmoothMeshElastic()
    {
        std::cout << "Elastic smoothing not (yet) implemented." << std::endl;
    }

};

}

#endif
