#ifndef DTCC_CONSTANTS_H
#define DTCC_CONSTANTS_H

class Constants
{
public:
  // Tolerance for geometric tests
  static constexpr double Epsilon = 1e-6;

  // Precision for output and printing
  static constexpr double Precision = 16;

  // Threshold for filtering points with small angles in building footprints
  static constexpr double FootprintAngleThreshold = 0.01;

  // Threshold for filtering outliers (clouds?) from point cloud
  static constexpr double PointCloudOutlierThreshold = 150.0;

  // Number of digits of precision used when writing files
  static constexpr double OutputPrecision = 3;
};

#endif