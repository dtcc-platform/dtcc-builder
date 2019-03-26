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

    //--- Public (run-time) parameters ---

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

protected:

    //--- Protected (compile-time) parameters ---

    // Tolerance for geometric tests
    static constexpr double Epsilon = 1e-6;

    // Threshold for filtering outliers (clouds?) from point cloud
    static constexpr double PointCloudOutlierThreshold = 150.0;

    // Classes with access to protected parameters
    friend class HeightMapGenerator;
    friend class CityModelGenerator;

};

std::ostream& operator<<(std::ostream& stream, const Parameters& parameters)
{
    stream << "Parameters:"
           << " XMin = "                << parameters.XMin
           << " YMin = "                << parameters.YMin
           << " XMax = "                << parameters.XMax
           << " YMax = "                << parameters.YMax
           << " DomainHeight = "        << parameters.DomainHeight
           << " HeightMapResolution = " << parameters.HeightMapResolution
           << " MeshResolution = "      << parameters.MeshResolution;
}

}

#endif
