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

  //--- Run-time parameters (parsed from file) ---

  // Origin
  double X0 = 0.0;
  double Y0 = 0.0;

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

  // Minimal building distance (merged if closer)
  double MinimalBuildingDistance = 0.5;

  //--- Compile-time parameters ---

  // Tolerance for geometric tests
  static constexpr double Epsilon = 1e-6;

  // Precision for output and printing
  static constexpr double Precision = 16;

  // Threshold for filtering duplicate points in building footprints
  static constexpr double FootprintDuplicateThreshold = 1.0;

  // Threshold for filtering outliers (clouds?) from point cloud
  static constexpr double PointCloudOutlierThreshold = 150.0;

};

std::ostream& operator<<(std::ostream& s, const Parameters& parameters)
{
  s << "Parameters:" << std::endl
    << "  X0                      = " << parameters.X0
    << std::endl
    << "  Y0                      = " << parameters.Y0
    << std::endl
    << "  XMin                    = " << parameters.XMin
    << std::endl
    << "  YMin                    = " << parameters.YMin
    << std::endl
    << "  XMax                    = " << parameters.XMax
    << std::endl
    << "  YMax                    = " << parameters.YMax
    << std::endl
    << "  DomainHeight            = " << parameters.DomainHeight
    << std::endl
    << "  HeightMapResolution     = " << parameters.HeightMapResolution
    << std::endl
    << "  MeshResolution          = " << parameters.MeshResolution
    << std::endl
    << "  MinimalBuildingDistance = " << parameters.MeshResolution;
}

}

#endif
