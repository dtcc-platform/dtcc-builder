// Copyright (C) 2020 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_POINT_CLOUD_PROCESSOR_H
#define DTCC_POINT_CLOUD_PROCESSOR_H

#include "Color.h"
#include "GeoRaster.h"
#include "PointCloud.h"

namespace DTCC
{

class PointCloudProcessor
{
public:
  // Georaster is assumed to have at least 3 bands where 1,2,3 are R,G,B
  // and in range 0-255
  static void ColorFromImage(PointCloud &pointCloud, const GeoRaster &raster )  
  {
      double r,g,b;
      pointCloud.Colors.clear();
      
      for (auto p: pointCloud.Points) {
        Point2D p2d = Point2D(p.x,p.y);
        try
        {
          r = raster(p2d,1);
          g = raster(p2d,2);
          b = raster(p2d,3);
          r /= 255;
          g /= 255;
          b /= 255;
        } catch (const std::runtime_error& error)
        {
          // point outside of raster
          r = 0.0;
          g = 0.0;
          b = 0.0;
        }
        
        pointCloud.Colors.push_back(Color(r,g,b));
      }

  }

  static PointCloud
  ClassificationFilter(const PointCloud &pointCloud,
                       const std::vector<int> &classifications)
  {
    PointCloud outCloud;
    size_t pointCount = pointCloud.Points.size();
    bool has_color = (pointCount == pointCloud.Colors.size());

    if (pointCount != pointCloud.Classification.size())
    {
      throw std::runtime_error("Point cloud isn't classified");
    }

    for (size_t i = 0; i < pointCount; i++)
    {
      for (auto c : classifications)
      {
        if (pointCloud.Classification[i] == c)
        {
          outCloud.Points.push_back(pointCloud.Points[i]);
          if (has_color)
          {
            outCloud.Colors.push_back(pointCloud.Colors[i]);
          }
          outCloud.Classification.push_back(pointCloud.Classification[i]);
          break;
        }
      }
    }
    return outCloud;
  }
};

} // namespace DTCC


#endif