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
  static PointCloud LoadPointCloud(std::string protobuf_string)
  {
    DTCC::PointCloud pc;
    pc.ParseFromString(protobuf_string);
    return LoadPointCloud(pc);
  }

  static PointCloud LoadPointCloud(DTCC::PointCloud pb_pc)
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

  static std::string ExportPointCloud(const PointCloud &pointCloud)
  {
    std::vector<DTCC::Vector3D> pb_pts;
    for (const auto &pt : pointCloud.Points)
    {
      pb_pts.push_back(DTCC::Point(pt.x, pt.y, pt.z));
    }
    std::vector<int> classifications;
    for (const auto c : pointCloud.Classifications)
      classifications.push_back(c);
    auto pb_pc = DTCC::CreatePointCloud(pb_pts, classifications);

    std::string pb_string;
    pb_pc.SerializeToString(&pb_string);

    return pb_string;
  }
};    
}
#endif