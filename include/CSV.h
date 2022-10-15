// Copyright (C) 2020-2021 Anders Logg, Vasilis Naserentin, Anton J Olsson
// Licensed under the MIT License

#ifndef DTCC_CSV_H
#define DTCC_CSV_H

#include "rapidcsv.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <vector>

#include "Mesh.h"
#include "Parameters.h"
#include "Vector.h"
#include "Simplex.h"
#include "PointCloud.h"
#include "Logging.h"

namespace DTCCBUILDER
{

class CSV
{
public:
  // Wrap for rapidcsv CSV container
  rapidcsv::Document Document;
  // Write 2D point set to CSV file
  static void Write(const std::vector<Vector2D> &Points, const std::string& fileName)
  {
    info("CSV: Writing 2D point set to file " + fileName);

    // Open file
    std::ofstream f(fileName.c_str());

    // Check file
    if (!f)
      throw std::runtime_error("Unable to write to file: " + fileName);

    // Set precision
    f << std::setprecision(Parameters::Precision);

    // Write points
    for (auto const &p : Points)
      Write(p, f);

    // Close file
    f.close();
  }

  // Write 2D mesh to CSV files (.csv will be appended)
  static void Write(const Mesh2D &mesh, const std::string& prefix)
  {
    info("CSV: Writing 2D mesh to file " + prefix);

    // Open files
    std::string fileNamePoints = prefix + "Points.csv";
    std::string fileNameTriangles = prefix + "Cells.csv";
    std::ofstream fp(fileNamePoints.c_str());
    std::ofstream ft(fileNameTriangles.c_str());

    // Check files
    if (!fp)
      throw std::runtime_error("Unable to write to file: " + fileNamePoints);
    if (!ft)
      throw std::runtime_error("Unable to write to file: " + fileNameTriangles);

    // Set precision
    fp << std::setprecision(Parameters::Precision);
    ft << std::setprecision(Parameters::Precision);

    // Write points
    for (auto const &p : mesh.Vertices)
      Write(Vector2D(p), fp);

    // Write triangles
    for (auto const &t : mesh.Cells)
      Write(t, ft);

    // Close files
    fp.close();
    ft.close();
  }

  // Write 3D mesh to CSV files
  static void Write(const Mesh3D &mesh, const std::string& prefix)
  {
    info("CSV: Writing 3D mesh to file " + prefix);

    // Open files
    std::string fileNamePoints = prefix + "Points.csv";
    std::string fileNameTriangles = prefix + "Cells.csv";
    std::ofstream fp(fileNamePoints.c_str());
    std::ofstream ft(fileNameTriangles.c_str());

    // Check files
    if (!fp)
      throw std::runtime_error("Unable to write to file: " + fileNamePoints);
    if (!ft)
      throw std::runtime_error("Unable to write to file: " + fileNameTriangles);

    // Write points
    for (auto const &p : mesh.Vertices)
      Write(Vector3D(p), fp);

    // Write triangles
    for (auto const &t : mesh.Cells)
      Write(t, ft);

    // Close files
    fp.close();
    ft.close();
  }

  /// Write PointCloud to CSV file.
  /// \param pointCloud PointCloud to read from
  /// \param fileName Name of file to write to
  static void Write(const PointCloud &pointCloud, const std::string &fileName)
  {
    std::cout << "CSV: "
              << "Writing pointcloud to file " << fileName << std::endl;

    // Open file
    std::ofstream f(fileName.c_str());

    // Check file
    if (!f)
      throw std::runtime_error("Unable to write to file: " + fileName);

    // Set precision
    f << std::setprecision(Parameters::Precision);

    // Write points
    if (pointCloud.Points.size() == pointCloud.Colors.size())
    {
      for (size_t i = 0; i < pointCloud.Points.size(); i++)
      {
        f << pointCloud.Points[i].x << "," << pointCloud.Points[i].y << ","
          << pointCloud.Points[i].z << "," << pointCloud.Colors[i].R << ","
          << pointCloud.Colors[i].G << "," << pointCloud.Colors[i].B << ","
          << pointCloud.Classifications[i] << std::endl;
      }
    }
    else
    {
      for (auto const &p : pointCloud.Points)
      {
        f << p.x << "," << p.y << "," << p.z << std::endl;
       }
    }

    // Close file
    f.close();
  }

  // Read CSV file and store it in rapidcsv
  void Read(const std::string &iFilename,
            rapidcsv::Document &doc,
            bool verbose = false)
  {
    // Open file and check
    try
    {
      doc.Load(iFilename);
    }
    catch (...)
    {
      throw std::runtime_error("Unable to read to file: " + iFilename);
    }

    if (verbose)
    {
      info("Read " + str(doc.GetColumnCount() + 1) +
           " columns"); // Columns seems to be 0-indexed counts
      info("Read " + str(doc.GetRowCount() + 1) + " rows");
    }
  }

  /// Read point cloud from CSV file. Note that points will be added
  /// to the given point cloud, enabling reading data from several
  /// CSV files into the same point cloud.
  /// \param pointCloud PointCloud to write to
  /// \param fileName CSV filename. First row is assumed to be labels and order
  /// of cols are assumed to be xyzrgbc (c = classification).
  static void Read(PointCloud &pointCloud, const std::string &fileName)
  {
    info(str("CSV: ") + str("Reading point cloud from file ") + fileName);
    Read_(pointCloud, fileName);
  }

  /// Read point cloud from CSV file only if they are within the BoundingBox.
  /// \param pointCloud PointCloud to write to
  /// \param fileName CSV filename. First row is assumed to be labels and order
  /// of cols are assumed to be xyzrgbc (c = classification).
  /// \param bbox BoundingBox2D of which the points must lie within
  static void Read(PointCloud &pointCloud,
                   const std::string &fileName,
                   const BoundingBox2D &bbox)
  {
    info("CSV: Reading point cloud from file: " + fileName + " bounded by " +
         str(bbox));
    Read_(pointCloud, fileName, &bbox);
  }

  /// Read point cloud from CSV file only if they have the defined
  /// classification
  /// \param pointCloud PointCloud to write to
  /// \param fileName CSV filename. First row is assumed to be labels and order
  /// of cols are assumed to be xyzrgbc (c = classification).
  /// \param classifications Vector of allowed classifications
  static void Read(PointCloud &pointCloud,
                   const std::string &fileName,
                   const std::vector<int> &classifications)
  {
    info("CSV: Reading point cloud from file: " + fileName +
         " with classifications " + GetClassString(classifications));
    Read_(pointCloud, fileName, nullptr, &classifications);
  }

  /// Read point cloud from CSV file only if they have the defined
  /// classification and are within the BoundingBox.
  /// \param pointCloud PointCloud to write to
  /// \param fileName CSV filename. First row is assumed to be labels and order
  /// of cols are assumed to be xyzrgbc (c = classification).
  /// \param classifications Vector of allowed classifications \param bbox
  /// BoundingBox2D of which the points must lie within
  static void Read(PointCloud &pointCloud,
                   const std::string &fileName,
                   const std::vector<int> &classifications,
                   const BoundingBox2D &bbox)
  {
    info("CSV: Reading point cloud from file: " + fileName + " bounded by " +
         str(bbox) + ", with classifications" +
         GetClassString(classifications));
    Read_(pointCloud, fileName, &bbox, &classifications);
  }

private:
  /// Construct string of allowed classifications.
  /// \param classifications Allowed classifications
  /// \return String of allowed classifications
  static std::string GetClassString(const std::vector<int> &classifications)
  {
    std::string classString;
    for (int i = 0; i < classifications.size(); ++i)
    {
      classString += std::to_string(classifications[i]) +
                     (i < classifications.size() - 1 ? ", " : "");
    }
    return classString;
  }

  /// Check if vector's x and y values are inside bounding box.
  /// \param bbox The 2D bounding box
  /// \param p The vector
  /// \return Whether the vector is inside bounding box
  static bool InsideBBox(const BoundingBox2D &bbox, const Vector3D &p)
  {
    return p.x >= bbox.P.x && p.x <= bbox.Q.x && p.y >= bbox.P.y &&
           p.y <= bbox.Q.y;
  }

  /// Check if point type/classification belongs to allowed classifications.
  /// \param classifications The allowed classifications
  /// \param type The point's type/classification
  /// \return Whether the point type/classification belongs to allowed
  /// classifications
  static bool HasRightClassification(const std::vector<int> &classifications,
                                     int type)
  {
    return std::find(classifications.begin(), classifications.end(), type) !=
           classifications.end();
  }

  /// Read points into PointCloud, possibly depending on certain constraints.
  /// \param pointCloud PointCloud to write to
  /// \param fileName CSV filename. First row is assumed to be labels and
  /// order of cols are assumed to be xyzrgbc (c = classification). \param
  /// bboxP BoundingBox2D of which the points must lie within \param classP
  /// The allowed point classifications
  static void Read_(PointCloud &pointCloud,
                    const std::string &fileName,
                    const BoundingBox2D *bboxP = nullptr,
                    const std::vector<int> *classP = nullptr)
  {
    std::regex numRe(".*\\d*");
    rapidcsv::Document doc(fileName, rapidcsv::LabelParams(-1, -1));

    for (int i = 1; i < doc.GetRowCount(); ++i)
    {
      Vector3D point;
      point.x = doc.GetCell<double>(0, i);
      point.y = doc.GetCell<double>(1, i);
      point.z = doc.GetCell<double>(2, i);

      int type = doc.GetCell<int>(6, i);

      // If point doesn't match some supplied constraint, continue iteration
      if (bboxP && !InsideBBox(*bboxP, point) ||
          classP && !HasRightClassification(*classP, type))
        continue;

      pointCloud.Points.push_back(point);
      pointCloud.Classifications.push_back(type);

      Color color;
      color.R = doc.GetCell<double>(3, i);
      color.G = doc.GetCell<double>(4, i);
      color.B = doc.GetCell<double>(5, i);
      pointCloud.Colors.push_back(color);
    }
  }

  // Write 2D point to file
  static void Write(const Vector2D &p, std::ofstream &f)
  {
    f << p.x << "," << p.y << std::endl;
  }

  // Write 3D point to file
  static void Write(const Vector3D &p, std::ofstream &f)
  {
    f << p.x << "," << p.y << "," << p.z << std::endl;
  }

  // Write 2D simplex to file
  static void Write(const Simplex2D &t, std::ofstream &f)
  {
    f << t.v0 << "," << t.v1 << "," << t.v2 << std::endl;
  }
  // Write 3D simplex to file
  static void Write(const Simplex3D &t, std::ofstream &f)
  {
    f << t.v0 << "," << t.v1 << "," << t.v2 << "," << t.v3 << std::endl;
  }
};

} // namespace DTCCBUILDER

#endif
