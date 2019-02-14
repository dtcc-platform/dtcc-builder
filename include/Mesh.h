// Mesh classes for 2D and 3D.
// Copyright (C) 2018 Anders Logg.

#ifndef VC_MESH_H
#define VC_MESH_H

#include <vector>
#include "Point.h"
#include "Simplex.h"

namespace VirtualCity
{

class Mesh2D
{
public:

    // List of points (vertices)
    std::vector<Point2D> Points;

    // List of cells (triangles)
    std::vector<Simplex2D> Cells;

    // List of domain markers
    std::vector<size_t> DomainMarkers;

    // Compute cell midpoint
    Point2D MidPoint(const Simplex2D& Cell) const
    {
        Point2D c;
        c += Points[Cell.v0];
        c += Points[Cell.v1];
        c += Points[Cell.v2];
        c /= 3.0;
        return c;
    }

};

class Mesh3D
{
public:

    // List of points (vertices)
    std::vector<Point3D> Points;

    // List of tetrahedra
    std::vector<Simplex3D> Cells;

    // List of domain markers
    std::vector<size_t> DomainMarkers;

    // Compute cell midpoint
    Point3D MidPoint(const Simplex3D& Cell) const
    {
        Point3D c;
        c += Points[Cell.v0];
        c += Points[Cell.v1];
        c += Points[Cell.v2];
        c += Points[Cell.v3];
        c /= 4.0;
        return c;
    }

};

std::ostream& operator<<(std::ostream& stream, const Mesh2D& m)
{
    stream << "2D mesh with " << m.Points.size() << " points and "
           << m.Cells.size() << " cells (triangles)";
}

std::ostream& operator<<(std::ostream& stream, const Mesh3D& m)
{
    stream << "3D mesh with " << m.Points.size() << " points and "
           << m.Cells.size() << " cells (tetrahedra)";
}

}

#endif
