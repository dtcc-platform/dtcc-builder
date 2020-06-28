// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_LAS_H
#define DTCC_LAS_H

#include <fstream>
#include <iostream>
#include <liblas/liblas.hpp>
#include <string>

#include "Color.h"
#include "Vector.h"
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
      const Vector3D p(_p.GetX(), _p.GetY(), _p.GetZ());
      
      liblas::Color const& color = _p.GetColor();
      const Color c(color.GetRed()/65535.0,color.GetGreen()/65535.0,color.GetBlue()/65535.0);

      // Update bounding box dimensions
      if (pointCloud.Points.size() == 0)
      {
        pointCloud.BoundingBox.P.x = p.x;
        pointCloud.BoundingBox.P.y = p.y;
        pointCloud.BoundingBox.Q.x = p.x;
        pointCloud.BoundingBox.Q.y = p.y;
      }
      else
      {
        pointCloud.BoundingBox.P.x = std::min(p.x, pointCloud.BoundingBox.P.x);
        pointCloud.BoundingBox.P.y = std::min(p.y, pointCloud.BoundingBox.P.y);
        pointCloud.BoundingBox.Q.x = std::max(p.x, pointCloud.BoundingBox.Q.x);
        pointCloud.BoundingBox.Q.y = std::max(p.y, pointCloud.BoundingBox.Q.y);
      }

      // Add point to point cloud
      pointCloud.Points.push_back(p);
      pointCloud.Colors.push_back(c);
    }
  }
};

} // namespace DTCC

#endif
