// CSV I/O
// Anders Logg 2018, Vasilis Naserentin 2019

#ifndef DTCC_CSV_H
#define DTCC_CSV_H

#include <fstream>
#include <iostream>
#include <vector>
#include "rapidcsv.h"

#include "Mesh.h"
#include "Parameters.h"
#include "Point.h"
#include "Simplex.h"

namespace DTCC
{

class CSV
{
public:
  // Wrap for rapidcsv CSV container
  rapidcsv::Document Document;
  // Write 2D point set to CSV file
  static void Write(const std::vector<Point2D> &Points, std::string fileName)
  {
    std::cout << "CSV: "
              << "Writing 2D point set to file " << fileName << std::endl;

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
  static void Write(const Mesh2D &mesh, std::string prefix)
  {
    std::cout << "CSV: "
              << "Writing 2D mesh to file " << prefix << std::endl;

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
      Write(p, fp);

    // Write triangles
    for (auto const &t : mesh.Cells)
      Write(t, ft);

    // Close files
    fp.close();
    ft.close();
  }

  // Write 3D mesh to CSV files
  static void Write(const Mesh3D &mesh, std::string prefix)
  {
    std::cout << "CSV: "
              << "Writing 3D mesh to file " << prefix << std::endl;

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
      Write(p, fp);

    // Write triangles
    for (auto const &t : mesh.Cells)
      Write(t, ft);

    // Close files
    fp.close();
    ft.close();
  }
  // Read CSV file and store it in rapidcsv
  void Read(std::string iFilename, bool verbose = false)
  {
    // Open file and check
    try
    {
      Document.Load(iFilename);
    }
    catch (...)
    {
      throw std::runtime_error("Unable to read to file: " + iFilename);
    }

    if (verbose)
    {
      std::cout << "Read " << Document.GetColumnCount() + 1
                << " columns" // seems to be 0-indexed counts
                << std::endl;
      std::cout << "Read " << Document.GetRowCount() + 1 << " rows"
                << std::endl;
    }
  }

private:
  // Write 2D point to file
  static void Write(const Point2D &p, std::ofstream &f)
  {
    f << p.x << "," << p.y << std::endl;
  }

  // Write 3D point to file
  static void Write(const Point3D &p, std::ofstream &f)
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

} // namespace DTCC

#endif
