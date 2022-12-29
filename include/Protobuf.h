// Copyright (C) 2022 Dag Wästberg
// Licensed under the MIT License
//
//  Modified by Anders Logg 2022

#ifndef DTCC_PROTOBUF
#define DTCC_PROTOBUF

#include "protobuf/dtcc.pb.h"
#include "protobuf/include/PointCloudMethods.h"
#include "protobuf/include/PolygonMethods.h"

#include "Logging.h"
#include "Point.h"
#include "PointCloud.h"
#include "PointCloudProcessor.h"
#include "Surface.h"

#include "datamodel/Building.h"
#include "datamodel/CityModel.h"

namespace DTCC_BUILDER
{

class Protobuf
{
public:
  // Developer note: Started but not completed migration to read/write functions
  // (issue #52) Developer note: DTCC Builder object / DTCC Protobuf _object

  // Read CityModel from Protobuf string
  static void read(CityModel &city_model, const std::string &pb)
  {
    DTCC::CityModel _city_model{};
    _city_model.ParseFromString(pb);
    read(city_model, _city_model);
  }

  // Read CityModel from Protobuf object
  static void read(CityModel &city_model, const DTCC::CityModel &_city_model)
  {
    for (const auto &_building : _city_model.buildings())
    {
      Building building{};
      building.UUID = _building.uuid();
      for (const auto &_vertex : _building.footprint().shell().vertices())
        building.Footprint.Vertices.push_back(
            Point2D(_vertex.x(), _vertex.y()));
      city_model.Buildings.push_back(building);
    }
  }

  // FIXME: Rename to void write()
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
      pb_bld.set_groundheight(building.GroundHeight);

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

    pb_cm.mutable_georeference()->set_x0(cityModel.Origin.x);
    pb_cm.mutable_georeference()->set_y0(cityModel.Origin.y);

    std::string pb_string;
    pb_cm.SerializeToString(&pb_string);
    // info("CityModel Serialized: " + str(pb_string.size()));

    return pb_string;
  }

  // FIXME: Rename to void read()
  static PointCloud LoadPointCloud(std::string protobuf_string)
  {
    DTCC::PointCloud pc;
    pc.ParseFromString(protobuf_string);
    return LoadPointCloud(pc);
  }

  // FIXME: Rename to void read()
  static PointCloud LoadPointCloud(DTCC::PointCloud pb_pc)
  {
    PointCloud pc;

    size_t num_pts = pb_pc.points_size();

    for (const auto &pts : pb_pc.points())
    {
      pc.Points.push_back(Point3D(pts.x(), pts.y(), pts.z()));
    }
    pc.CalculateBoundingBox();

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

    if (pb_pc.returnnumber_size() > 0)
    {
      if (pb_pc.returnnumber_size() != num_pts ||
          pb_pc.numreturns_size() != num_pts)
        error("number of scan flags not equal to number of points!");
      for (size_t i = 0; i < num_pts; i++)
      {
        pc.ScanFlags.push_back(PointCloudProcessor::packScanFlag(
            pb_pc.returnnumber(i), pb_pc.numreturns(i)));
      }
    }

    return pc;
  }

  // FIXME: Rename to void write()
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

    pb_pc.mutable_georeference()->set_x0(pointCloud.Origin.x);
    pb_pc.mutable_georeference()->set_y0(pointCloud.Origin.y);

    std::string pb_string;
    pb_pc.SerializeToString(&pb_string);

    return pb_string;
  }

  // Read Surface3D from Protobuf string
  static void read(Surface3D &surface, const std::string &pb)
  {
    error("Reading Surface3D from Protobuf not implemented");
  }

  // Read Surface3D from Protobuf object
  static void read(Surface3D &surface, const DTCC::Surface3D &_surface)
  {
    error("Reading Surface3D from Protobuf not implemented");
  }

  // Write Surface3D to Protobuf string
  static void write(const Surface3D &surface, std::string &pb)
  {
    DTCC::Surface3D _surface{};

    // Write vertices
    for (const auto &vertex : surface.Vertices)
    {
      auto _vertex = _surface.add_vertices();
      _vertex->set_x(vertex.x);
      _vertex->set_y(vertex.y);
      _vertex->set_z(vertex.z);
    }

    // Write faces
    for (const auto &face : surface.Faces)
    {
      auto _face = _surface.add_faces();
      _face->set_v0(face.v0);
      _face->set_v1(face.v1);
      _face->set_v2(face.v2);
    }

    // Write normals
    for (const auto &normal : surface.Normals)
    {
      auto _normal = _surface.add_normals();
      _normal->set_x(normal.x);
      _normal->set_y(normal.y);
      _normal->set_z(normal.z);
    }

    // Write to Protobuf string
    _surface.SerializeToString(&pb);
  }
};
} // namespace DTCC_BUILDER
#endif
