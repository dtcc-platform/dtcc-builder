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

    // Radius of computational domain relative to radius of city model
    double DomainRadius = 2.0;

    // Maximum mesh size used for mesh generation [m]
    double MeshSize = 1.0;

};

}

#endif
