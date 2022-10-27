// Copyright (C) 2022 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_PROTOBUF_CONVERTER
#define DTCC_PROTOBUF_CONVERTER

#include "protobuf/dtcc.pb.h"
#include "protobuf/include/PointCloudMethods.h"
#include "protobuf/include/PolygonMethods.h"

#include "Logging.h"
#include "Point.h"
#include "PointCloud.h"
#include "PointCloudProcessor.h"
namespace DTCC_BUILDER
{

class ProtobufConverter
{
public:
  static PointCloud ConvertPointCloud(DTCC::PointCloud pb_pc)
  {
    PointCloud pc;

    size_t num_pts = pb_pc.points_size();

    for (const auto &pts : pb_pc.points())
    {
      pc.Points.push_back(Point3D(pts.x(), pts.y(), pts.z()));
    }

    if (pb_pc.classification_size() > 0)
    {
      if (pb_pc.classification_size() != num_pts)
        error("number of classification not equal to number of points!");
      for (const auto &cl : pb_pc.classification())
      {
        pc.Classifications.push_back(cl);
      }
    }

    if (pb_pc.intensity_size() > 0)
    {
      for (const auto &i : pb_pc.intensity())
      {
        pc.Intensities.push_back(i);
      }
    }

    if (pb_pc.returnnumber_size()>0)
    {
      if (pb_pc.returnnumber_size() != num_pts || 
          pb_pc.numreturns_size() != num_pts)
        error("number of scan flags not equal to number of points!");
      for (size_t i = 0; i < num_pts; i++)
      {
        pc.ScanFlags.push_back(PointCloudProcessor::packScanFlag(
          pb_pc.returnnumber(i) ,pb_pc.numreturns(i)));
      }
      
    }
    
    return pc;
  }
};    
}
#endif