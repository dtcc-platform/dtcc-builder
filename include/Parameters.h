// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_PARAMETERS_H
#define DTCC_PARAMETERS_H

#include <string>
#include "Logging.h"

namespace DTCC
{

  class Parameters : public Printable
  {
  public:
    //--- Run-time parameters (parsed from file) ---

    // Directory for input/output
    std::string DataDirectory;

    // Origin
    double X0 = 0.0;
    double Y0 = 0.0;

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
    double MinimalBuildingDistance = 1.0;

    // Minimal vertex distance (merged if closer)
    double MinimalVertexDistance = 1.0;

    // Height of computational domain
    double DomainHeight = 100.0;

    // Maximum mesh size used for mesh generation [m]
    double MeshResolution = 10.0;

    // Keep ground flat (ignore height map)
    bool FlatGround = false;

    // Number of smoothing iterations for extra ground mesh
    int GroundSmoothing = 5;

    // Write extra data for debugging
    bool Debug = false;

    //--- Compile-time parameters ---

    // Tolerance for geometric tests
    static constexpr double Epsilon = 1e-6;

    // Precision for output and printing
    static constexpr double Precision = 16;

    // Threshold for filtering points with small angles in building footprints
    static constexpr double FootprintAngleThreshold = 0.01;

    // Threshold for filtering outliers (clouds?) from point cloud
    static constexpr double PointCloudOutlierThreshold = 150.0;

    /// Pretty-print
    std::string __str__() const override
    {
      return str("Parameters:") +
             "\n  DataDirectory           = " + DataDirectory +
             "\n  X0                      = " + str(X0) +
             "\n  Y0                      = " + str(Y0) +
             "\n  XMin                    = " + str(XMin) +
             "\n  YMin                    = " + str(YMin) +
             "\n  XMax                    = " + str(XMax) +
             "\n  YMax                    = " + str(YMax) +
             "\n  AutoDomain              = " + str(AutoDomain) +
             "\n  HeightMapResolution     = " + str(HeightMapResolution) +
             "\n  MinimalBuildingDistance = " + str(MinimalBuildingDistance) +
             "\n  DomainHeight            = " + str(DomainHeight) +
             "\n  MeshResolution          = " + str(MeshResolution) +
             "\n  FlatGround              = " + str(FlatGround);
    }
  };

} // namespace DTCC

#endif
