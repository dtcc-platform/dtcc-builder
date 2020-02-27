// LAS I/O
// Anders Logg 2019

#ifndef DTCC_SHP_H
#define DTCC_SHP_H

#include <fstream>
#include <iostream>
#include <liblas/liblas.hpp>
#include <string>

#include "Point.h"
#include "PointCloud.h"

namespace DTCC
{

class LAS
{
public:
  // Read point cloud from LAS file. Note that points will be added
  // to the given point cloud, enabling reading data from several
  // LAS files into the same point cloud.
  static void Read(PointCloud &pointCloud, std::string fileName)
  {
    std::cout << "LAS: "
              << "Reading point cloud from file " << fileName << std::endl;

    // Open file
    std::ifstream f;
    f.open(fileName, std::ios::in | std::ios::binary);

    // Create reader
    liblas::ReaderFactory factory;
    liblas::Reader reader = factory.CreateWithStream(f);

    // Read header
    liblas::Header const &header = reader.GetHeader();
    const bool isCompressed = header.Compressed();
    const std::string signature = header.GetFileSignature();
    const size_t numPoints = header.GetPointRecordsCount();
    if (isCompressed)
      std::cout << "LAS: Compressed" << std::endl;
    else
      std::cout << "LAS: Uncompressed" << std::endl;
    std::cout << "LAS: " << signature << std::endl;
    std::cout << "LAS: " << numPoints << " points" << std::endl;

    // Iterate over points
    while (reader.ReadNextPoint())
    {
      // Get point
      liblas::Point const &_p = reader.GetPoint();
      const Point3D p(_p.GetX(), _p.GetY(), _p.GetZ());

      // Update bounding box dimensions
      if (pointCloud.Points.size() == 0)
      {
        pointCloud.XMin = p.x;
        pointCloud.YMin = p.y;
        pointCloud.XMax = p.x;
        pointCloud.YMax = p.y;
      }
      else
      {
        pointCloud.XMin = std::min(p.x, pointCloud.XMin);
        pointCloud.YMin = std::min(p.y, pointCloud.YMin);
        pointCloud.XMax = std::max(p.x, pointCloud.XMax);
        pointCloud.YMax = std::max(p.y, pointCloud.YMax);
      }

      // Add point to point cloud
      pointCloud.Points.push_back(p);
    }
  }
};

} // namespace DTCC

#endif
