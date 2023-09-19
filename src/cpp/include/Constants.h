#ifndef DTCC_CONSTANTS_H
#define DTCC_CONSTANTS_H

class Constants
{
public:
  // Tolerance for geometric tests
  static constexpr double epsilon = 1e-6;

  // precision for output and printing
  static constexpr double precision = 16;

  // Threshold for filtering points with small angles in building footprints
  static constexpr double footprint_angle_threshold = 0.01;

  // Threshold for filtering outliers (clouds?) from point cloud
  static constexpr double point_cloud_outlier_threshold = 150.0;

  // Number of digits of precision used when writing files
  static constexpr double output_precision = 3;
};

#endif