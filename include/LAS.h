// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_LAS_H
#define DTCC_LAS_H

#include <fstream>
#include <iostream>
#include <liblas/liblas.hpp>
#include <string>

#include "BoundingBox.h"
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
    std::vector<liblas::FilterPtr> filters;
    _Read(pointCloud,fileName,filters);
  }

  // Read point cloud from LAS file only if they are within the BoundingBox
  static void Read(PointCloud &pointCloud, std::string fileName, BoundingBox2D bbox )
  {
    std::cout << "LAS: Reading point cloud from file: " << fileName 
              << " bounded by " << str(bbox) << std::endl;

    liblas::Bounds<double> bounds;
    std::vector<liblas::FilterPtr> filters;
    bounds = liblas::Bounds<double>(bbox.P.x,bbox.P.y,bbox.Q.x,bbox.Q.y);
    liblas::FilterPtr bounds_filter = liblas::FilterPtr(new liblas::BoundsFilter(bounds));
    bounds_filter->SetType(liblas::FilterI::eInclusion);
    filters.push_back(bounds_filter);

    _Read(pointCloud, fileName, filters);
  }

  // Read point cloud from LAS file only if they have the defined classification
  static void Read(PointCloud &pointCloud, std::string fileName, const std::vector<int> &classifications) 
  {
    std::vector<liblas::Classification> classes;
    for (int c: classifications ) 
    {
      classes.push_back(liblas::Classification(c));
    }
    std::vector<liblas::FilterPtr> filters;
    liblas::FilterPtr class_filter =
        liblas::FilterPtr(new liblas::ClassificationFilter(classes));
    class_filter->SetType(liblas::FilterI::eInclusion);
    filters.push_back(class_filter);

    _Read(pointCloud, fileName, filters);
  }

  // Read point cloud from LAS file only if they have the defined classification and are within the BoundingBox
  static void Read(PointCloud &pointCloud, std::string fileName, const std::vector<int> &classifications,BoundingBox2D bbox) 
  {
    std::cout << "LAS: Reading point cloud from file: " << fileName 
              << " bounded by " << str(bbox) << std::endl;
    std::vector<liblas::Classification> classes;
    liblas::Bounds<double> bounds;
    std::vector<liblas::FilterPtr> filters;


    for (int c: classifications ) 
    {
      classes.push_back(liblas::Classification(c));
    }
    liblas::FilterPtr class_filter =
    liblas::FilterPtr(new liblas::ClassificationFilter(classes));
    class_filter->SetType(liblas::FilterI::eInclusion);
    filters.push_back(class_filter);

    bounds = liblas::Bounds<double>(bbox.P.x,bbox.P.y,bbox.Q.x,bbox.Q.y);
    liblas::FilterPtr bounds_filter = liblas::FilterPtr(new liblas::BoundsFilter(bounds));
    bounds_filter->SetType(liblas::FilterI::eInclusion);
    filters.push_back(bounds_filter);

    _Read(pointCloud, fileName, filters);
  }



private:
  static void _Read(PointCloud &pointCloud, std::string fileName,std::vector<liblas::FilterPtr> filters) 
  {
    // Open file
    std::ifstream f;
    f.open(fileName, std::ios::in | std::ios::binary);
    

    // Create reader
    liblas::ReaderFactory factory;
    liblas::Reader reader = factory.CreateWithStream(f);

    if (filters.size() > 0) {
      reader.SetFilters(filters);
    }
    // Read header
    liblas::Header const &header = reader.GetHeader();
    const bool isCompressed = header.Compressed();
    const std::string signature = header.GetFileSignature();
    const size_t numPoints = header.GetPointRecordsCount();
    size_t readPoints = 0;
    if (isCompressed)
      std::cout << "LAS: Compressed" << std::endl;
    else
      std::cout << "LAS: Uncompressed" << std::endl;
    std::cout << "LAS: " << signature << std::endl;
    std::cout << "LAS: contains " << numPoints << " points" << std::endl;

    // Iterate over points
    while (reader.ReadNextPoint())
    {
      // Get point
      liblas::Point const &_p = reader.GetPoint();
      const Vector3D p(_p.GetX(), _p.GetY(), _p.GetZ());

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
      readPoints++;
    }
    std::cout << "LAS: read " << readPoints << " points" << std::endl;

  }

};

} // namespace DTCC

#endif
