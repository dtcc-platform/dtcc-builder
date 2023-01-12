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

#include "datamodel/Building.h"
#include "datamodel/CityModel.h"

namespace DTCC_BUILDER
{

class ProtobufConverter
{
public:
  static CityModel LoadCityModel(std::string protobuf_string)
  {
    DTCC::CityModel cm;
    cm.ParseFromString(protobuf_string);
    return LoadCityModel(cm);
  }

  static CityModel LoadCityModel(DTCC::CityModel pb_cm)
  {
    CityModel cm;
    auto pb_buildings = pb_cm.buildings();
    for (const auto &pb_building : pb_buildings)
    {
      Building b;
      b.UUID = pb_building.uuid();
      auto pb_footprint = pb_building.footprint().shell().vertices();
      for (const auto &v : pb_footprint)
      {
        b.Footprint.Vertices.push_back(Point2D(v.x(), v.y()));
      }

      cm.Buildings.push_back(b);
    }

    return cm;
  }

  static std::string ExportCityModel(const CityModel &cityModel)
  {
    DTCC::CityModel pb_cm;

    std::vector<DTCC::Building> pb_buildings;
    info("Serializing " + str(cityModel.Buildings.size()) + " buildings");
    for (const auto &building : cityModel.Buildings)
    {
      DTCC::Building pb_bld;
      pb_bld.set_height(building.Height);
      pb_bld.set_uuid(building.UUID);
      pb_bld.set_error(building.error);

      std::vector<DTCC::Vector2D> pb_verts;
      for (const auto &vert : building.Footprint.Vertices)
      {
        pb_verts.push_back(DTCC::Vertex(vert.x, vert.y));
      }

      pb_bld.mutable_footprint()->CopyFrom(CreatePolygon(pb_verts));
      pb_buildings.push_back(pb_bld);
    }

    google::protobuf::RepeatedPtrField<DTCC::Building> building_data(
        pb_buildings.begin(), pb_buildings.end());

    pb_cm.mutable_buildings()->Swap(&building_data);

    std::string pb_string;
    pb_cm.SerializeToString(&pb_string);
    info("CityModel Serialized: " + str(pb_string.size()));

    return pb_string;
  }

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