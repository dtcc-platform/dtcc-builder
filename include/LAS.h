// Copyright (C) 2020 Anders Logg, Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_LAS_H
#define DTCC_LAS_H

#include <fstream>
#include <iostream>
#include <liblas/liblas.hpp>
#include <string>
#include <utility>

#include "BoundingBox.h"
#include "Color.h"
#include "CommandLine.h"
#include "Filesystem.h"
#include "PointCloud.h"
#include "Vector.h"

namespace DTCC
{

class LAS
{
public:
  // Read point cloud from LAS file. Note that points will be added
  // to the given point cloud, enabling reading data from several
  // LAS files into the same point cloud.
  // if extraData is true then information about intensity and scan returns will
  // also be read (if available)
  static void Read(PointCloud &pointCloud,
                   const std::string &fileName,
                   bool extraData = false)
  {
    Info(str("LAS: ") + str("Reading point cloud from file ") + fileName);
    std::vector<liblas::FilterPtr> filters;
    _Read(pointCloud, fileName, filters, extraData);
  }

  /// Read point cloud from LAS file and filter by bounding box.
  static void Read(PointCloud &pointCloud,
                   const std::string &fileName,
                   const BoundingBox2D &bbox,
                   bool extraData = false)
  {
    Info("LAS: Reading point cloud from file: " + fileName + " bounded by " +
         str(bbox));

    liblas::Bounds<double> bounds;
    std::vector<liblas::FilterPtr> filters;
    filters.push_back(MakeBoundsFilter(bbox));
    _Read(pointCloud, fileName, filters, extraData);
  }

  /// Read point cloud from LAS file and filter by classification.
  static void Read(PointCloud &pointCloud,
                   const std::string &fileName,
                   const std::vector<int> &classifications,
                   bool extraData = false)
  {
    std::vector<liblas::FilterPtr> filters;
    filters.push_back(MakeClassFilter(classifications));
    _Read(pointCloud, fileName, filters, extraData);
  }

  /// Read point cloud from LAS file and filter by classification and bounding
  /// box.
  static void Read(PointCloud &pointCloud,
                   const std::string &fileName,
                   const std::vector<int> &classifications,
                   const BoundingBox2D &bbox,
                   bool extraData = false)
  {

    Info("LAS: Reading point cloud from file: " + fileName + " bounded by " +
         str(bbox));
    std::vector<liblas::FilterPtr> filters;
    filters.push_back(MakeClassFilter(classifications));
    filters.push_back(MakeBoundsFilter(bbox));
    _Read(pointCloud, fileName, filters, extraData);
  }

  /// Read point cloud from directory of LAS files (.las and .laz).
  static void ReadDirectory(PointCloud &pointCloud,
                            const std::string &directoryName,
                            bool extraData = false)
  {
    std::string dir{directoryName};
    if (!Utils::EndsWith(dir, "/"))
      dir += "/";
    Info("Reading all .las and .laz files from directory " + dir + "...");
    for (auto const &f : Filesystem::ListDirectory(dir))
      if (Utils::EndsWith(f, ".las") || Utils::EndsWith(f, ".laz"))
        Read(pointCloud, f, extraData);
  }

  /// Read point cloud from directory of LAS files (.las and .laz) and filter by
  /// bounding box
  static void ReadDirectory(PointCloud &pointCloud,
                            const std::string &directoryName,
                            const BoundingBox2D &bbox,
                            bool extraData = false)
  {
    std::string dir{directoryName};
    if (!Utils::EndsWith(dir, "/"))
      dir += "/";
    Info("Reading all .las and .laz files from directory " + dir + "...");
    for (auto const &f : Filesystem::ListDirectory(dir))
      if (Utils::EndsWith(f, ".las") || Utils::EndsWith(f, ".laz"))
        Read(pointCloud, f, bbox, extraData);
  }

  static void Bounds(BoundingBox2D &bb, const std::string &fileName)
  {
    // Open file
    std::ifstream f;
    f.open(fileName, std::ios::in | std::ios::binary);

    // Create reader
    liblas::ReaderFactory factory;
    liblas::Reader reader = factory.CreateWithStream(f);
    // Read header
    liblas::Header const &header = reader.GetHeader();

    auto bounds = header.GetExtent();
    bb.P.x = bounds.minx();
    bb.Q.x = bounds.maxx();
    bb.P.y = bounds.miny();
    bb.Q.y = bounds.maxy();
  }

  static void BoundsDirectory(BoundingBox2D &bb,
                              const std::string &directoryName)
  {
    std::string dir{directoryName};
    if (!Utils::EndsWith(dir, "/"))
      dir += "/";
    bool initBB = false;
    for (auto const &f : Filesystem::ListDirectory(dir))
    {
      if (Utils::EndsWith(f, ".las") || Utils::EndsWith(f, ".laz"))
      {
        auto fileBB = BoundingBox2D();
        LAS::Bounds(fileBB, f);
        if (!initBB)
        {
          bb.P = fileBB.P;
          bb.Q = fileBB.Q;
          initBB = true;
        }
        else
        {
          bb.Union(fileBB);
        }
      }
    }
  }

  /// Write point cloud to file
  static void Write(const PointCloud &pointCloud, const std::string &fileName)
  {
    Info("LAS: Writing " + str(pointCloud.Points.size()) + " points to " +
         fileName);

    std::ofstream ofs;
    ofs.open(fileName, std::ios::out | std::ios::binary);

    liblas::Header header;
    liblas::Writer writer(ofs, header);
    if (pointCloud.Points.size() == pointCloud.Colors.size())
    {
      for (size_t i = 0; i < pointCloud.Points.size(); i++)
      {
        liblas::Point point(&header);
        point.SetCoordinates(pointCloud.Points[i].x, pointCloud.Points[i].y,
                             pointCloud.Points[i].z);
        point.SetColor(liblas::Color(pointCloud.Colors[i].R * 65535,
                                     pointCloud.Colors[i].G * 65535,
                                     pointCloud.Colors[i].B * 65535));
        writer.WritePoint(point);
      }
    }
    else
    {
      for (auto const &p : pointCloud.Points)
      {
        liblas::Point point(&header);
        point.SetCoordinates(p.x, p.y, p.z);
        writer.WritePoint(point);
      }
    }
  }

private:
  static liblas::FilterPtr
  MakeClassFilter(const std::vector<int> &classifications)
  {
    std::vector<liblas::Classification> classes;
    for (int c : classifications)
    {
      classes.push_back(liblas::Classification(c));
    }
    liblas::FilterPtr class_filter =
        liblas::FilterPtr(new liblas::ClassificationFilter(classes));
    class_filter->SetType(liblas::FilterI::eInclusion);

    return class_filter;
  }

  static liblas::FilterPtr MakeBoundsFilter(const BoundingBox2D &bbox)
  {
    liblas::Bounds<double> bounds;
    bounds = liblas::Bounds<double>(bbox.P.x, bbox.P.y, bbox.Q.x, bbox.Q.y);
    liblas::FilterPtr bounds_filter =
        liblas::FilterPtr(new liblas::BoundsFilter(bounds));
    bounds_filter->SetType(liblas::FilterI::eInclusion);

    return bounds_filter;
  }

  static void _Read(PointCloud &pointCloud,
                    const std::string &fileName,
                    const std::vector<liblas::FilterPtr> &filters,
                    bool extraData)
  {
    // Open file
    std::ifstream f;
    f.open(fileName, std::ios::in | std::ios::binary);

    // Create reader
    liblas::ReaderFactory factory;
    liblas::Reader reader = factory.CreateWithStream(f);
    if (!filters.empty())
      reader.SetFilters(filters);

    // Read header
    liblas::Header const &header = reader.GetHeader();
    const bool isCompressed = header.Compressed();
    const std::string signature = header.GetFileSignature();
    const size_t numPoints = header.GetPointRecordsCount();
    const std::string compression =
        (isCompressed ? "Compressed" : "Uncompressed");
    Info("LAS: " + compression + ", " + signature + ", " + str(numPoints) +
         " points");

    // Iterate over points
    size_t readPoints = 0;
    while (reader.ReadNextPoint())
    {
      // Get point
      liblas::Point const &_p = reader.GetPoint();
      const Vector3D p(_p.GetX(), _p.GetY(), _p.GetZ());

      liblas::Color const &color = _p.GetColor();
      liblas::Classification const &cls = _p.GetClassification();
      uint16_t intensity = 0;
      uint8_t scanFlags = 0;
      if (extraData)
      {
        intensity = _p.GetIntensity();
        scanFlags = _p.GetScanFlags();
      }

      // Colors seem to be 16-bit in las spec.
      const Color c(color.GetRed() / 65535.0, color.GetGreen() / 65535.0,
                    color.GetBlue() / 65535.0);

      // Update bounding box dimensions
      if (pointCloud.Points.empty())
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
      pointCloud.Classifications.push_back(cls.GetClass());
      if (extraData)
      {
        pointCloud.Intensities.push_back(intensity);
        pointCloud.ScanFlags.push_back(scanFlags);
      }
      readPoints++;
    }

    Info("LAS: Found " + str(readPoints) + " points matching given filters");
  }
};
} // namespace DTCC

#endif
