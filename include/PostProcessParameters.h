// Post process parameters for dtccore.
// Anders Logg, Vasilis Naserentin 2020
// Licensed under the MIT License

#ifndef DTCC_POST_PROCESS_PARAMETERS_H
#define DTCC_POST_PROCESS_PARAMETERS_H

#include <string>
#include <vector>

namespace DTCC
{

class PostProcessVariables

{

public:
  std::string name = "";
  int sizex = 0;
  int sizey = 0;
  int sizez = 0;
};

class PostProcessParameters
{
public:
  //--- Run-time parameters (parsed from file) ---

  // Directory for input/output
  std::string DataDirectory;

  // Origin
  double X0 = 0.0;
  double Y0 = 0.0;
  // Variable
  std::vector<PostProcessVariables> Variables;
  std::string Variable;
  /*
  // Domain dimensions
  double XMin = 0.0;
  double YMin = 0.0;
  double XMax = 0.0;
  double YMax = 0.0;

  // Automatically determine domain size
  bool AutoDomain = true;

  // Height map resolution
  double HeightMapResolution = 1.0;

  // Minimal building distance (merged if closer)
  double MinimalBuildingDistance = 0.5;

  // Height of computational domain
  double DomainHeight = 100.0;

  // Maximum mesh size used for mesh generation [m]
  double MeshResolution = 10.0;

  // Keep ground flat (ignore height map)
  bool FlatGround = false;

  // Number of smoothing iterations for extra ground mesh
  int GroundSmoothing = 5;

  //--- Compile-time parameters ---

  // Tolerance for geometric tests
  static constexpr double Epsilon = 1e-6;

  // Precision for output and printing
  static constexpr double Precision = 16;

  // Threshold for filtering duplicate points in building footprints
  static constexpr double FootprintDuplicateThreshold = 1.0;

  // Threshold for filtering outliers (clouds?) from point cloud
  static constexpr double PointCloudOutlierThreshold = 150.0;
  */
};

std::ostream &operator<<(std::ostream &s, const PostProcessVariables &variables)
{
  s << "Post Process Variables:" << variables.name << std::endl;

  return s;
}

std::ostream &operator<<(std::ostream &s,
                         const std::vector<PostProcessVariables> &variables)
{
  s << "Post Process Variables:" << variables[0].name << std::endl;

  return s;
}

std::ostream &operator<<(std::ostream &s,
                         const PostProcessParameters &parameters)
{
  s << "Post Process Parameters:" << std::endl
    << "  X0                      = " << parameters.X0 << std::endl
    << "  Y0                      = " << parameters.Y0 << std::endl
    << "  Variable                = " << parameters.Variable << std::endl;

  return s;
}

} // namespace DTCC

#endif
