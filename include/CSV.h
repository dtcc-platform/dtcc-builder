// CSV I/O
// Anders Logg 2018

#ifndef CSV_H
#define CSV_H

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>

#include "Point.h"
#include "Simplex.h"
#include "Mesh.h"

namespace VirtualCity
{

class CSV
{
public:

    static const int PRECISION = 16;

    // Write 2D point set to CSV file (.csv will be appended)
    static void Write(const std::vector<Point2D>& Points,
                      std::string Prefix)
    {
        std::cout << "CSV: " << "Writing 2D point set to file "
                  << Prefix << std::endl;

        // Open file
        std::string FileName = Prefix + ".csv";
        std::ofstream f(FileName.c_str());

        // Check file
        if (!f)
            throw std::runtime_error("Unable to write to file: " + FileName);

        // Set precision
        f << std::setprecision(PRECISION);

        // Write points
        for (auto const & p : Points)
            Write(p, f);

        // Close file
        f.close();
    }

    // Write 2D mesh to CSV files (.csv will be appended)
    static void Write(const Mesh2D& Mesh,
                      std::string Prefix)
    {
        std::cout << "CSV: " << "Writing 2D mesh set to file "
                  << Prefix << std::endl;

        // Open files
        std::string FileNamePoints = Prefix + "Points.csv";
        std::string FileNameTriangles = Prefix + "Cells.csv";
        std::ofstream fp(FileNamePoints.c_str());
        std::ofstream ft(FileNameTriangles.c_str());

        // Check files
        if (!fp)
            throw std::runtime_error("Unable to write to file: " + FileNamePoints);
        if (!ft)
            throw std::runtime_error("Unable to write to file: " + FileNameTriangles);

        // Set precision
        fp << std::setprecision(PRECISION);
        ft << std::setprecision(PRECISION);

        // Write points
        for (auto const & p : Mesh.Points)
            Write(p, fp);

        // Write triangles
        for (auto const & t : Mesh.Cells)
            Write(t, ft);

        // Close files
        fp.close();
        ft.close();
    }

    // Write 3D mesh to CSV files (.csv will be appended)
    static void Write(const Mesh3D& Mesh,
                      std::string Prefix)
    {
        std::cout << "CSV: " << "Writing 3D mesh to file "
                  << Prefix << std::endl;

        // Open files
        std::string FileNamePoints = Prefix + "Points.csv";
        std::string FileNameTriangles = Prefix + "Cells.csv";
        std::ofstream fp(FileNamePoints.c_str());
        std::ofstream ft(FileNameTriangles.c_str());

        // Check files
        if (!fp)
            throw std::runtime_error("Unable to write to file: " + FileNamePoints);
        if (!ft)
            throw std::runtime_error("Unable to write to file: " + FileNameTriangles);

        // Write points
        for (auto const & p : Mesh.Points)
            Write(p, fp);

        // Write triangles
        for (auto const & t : Mesh.Cells)
            Write(t, ft);

        // Close files
        fp.close();
        ft.close();
    }

private:

    // Write 2D point to file
    static void Write(const Point2D& p, std::ofstream& f)
    {
        f << p.x << "," << p.y << std::endl;
    }

    // Write 3D point to file
    static void Write(const Point3D& p, std::ofstream& f)
    {
        f << p.x << "," << p.y << "," << p.z << std::endl;
    }

    // Write 2D simplex to file
    static void Write(const Simplex2D& t, std::ofstream& f)
    {
        f << t.v0 << "," << t.v1 << "," << t.v2 << std::endl;
    }

    // Write 3D simplex to file
    static void Write(const Simplex3D& t, std::ofstream& f)
    {
        f << t.v0 << "," << t.v1 << "," << t.v2 << "," << t.v3 << std::endl;
    }

};

}

#endif
