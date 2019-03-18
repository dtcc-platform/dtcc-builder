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

    // Domain dimensions
    double XMin = 0.0;
    double YMin = 0.0;
    double XMax = 100.0;
    double YMax = 100.0;

    // Height of computational domain
    double DomainHeight = 100.0;

    // Height map resolution
    double HeightMapResolution = 1.0;

    // Maximum mesh size used for mesh generation [m]
    double MeshResolution = 10.0;

};

}

#endif
