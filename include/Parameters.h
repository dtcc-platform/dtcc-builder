// Global parameters for VCCore.
// Anders Logg 2019

#ifndef VC_PARAMETERS_H
#define VC_PARAMETERS_H

#include <string>

namespace VirtualCity
{

class Parameters
{
public:

    // Center of computational domain
    double DomainCenterX = -3000.0;
    double DomainCenterY = -4000.0;

    // Radius of computational domain
    double DomainRadius = 100.0;

    // Height of computational domain
    double DomainHeight = 100.0;

    // Maximum mesh size used for mesh generation [m]
    double MeshSize = 10.0;

    // Stride length used for downsampling height map
    size_t HeightMapStride = 1;

};

}

#endif
