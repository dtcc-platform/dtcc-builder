// Representation of 2D height maps.
// Copyright (C) 2019 Anders Logg.

#ifndef VC_HEIGHT_MAP_H
#define VC_HEIGHT_MAP_H

#include <vector>

#include "Point.h"
#include "GeoReference.h"

namespace VirtualCity
{

class HeightMap
{
public:

    // Grid width
    size_t Width;

    // Grid height
    size_t Height;

    // Grid data (flattened array of (x, y) coordinates)
    std::vector<double> GridData;

    // Grid map (transform of grid data)
    GeoReference GridMap;

    // Create empty height map
    HeightMap() : Width(0), Height(0) {}

    // Return height (z) at world (UTM) point p
    double operator() (const Point2D& p) const
    {
        return (*this)(p.x, p.y);
    }

    // Return height (z) at world (UTM) coordinates (x, y)
    double operator() (double x, double y) const
    {
        // Transform world (UTM) coordinates to pixel coordinates
        const double Ex = GridMap.E * x;
        const double By = GridMap.B * y;
        const double BF = GridMap.B * GridMap.F;
        const double EC = GridMap.E * GridMap.C;
        const double Dx = GridMap.D * x;
        const double Ay = GridMap.A * y;
        const double DC = GridMap.D * GridMap.C;
        const double AF = GridMap.A * GridMap.F;
        const double AE = GridMap.A * GridMap.E;
        const double DB = GridMap.D * GridMap.B;
        const double det = AE - DB;
        double X = (Ex - By + BF - EC) / det;  // column
        double Y = (-Dx + Ay + DC - AF) / det; // row

        // Find the square containing the coordinate
        const std::size_t X0 = std::lround(X);
        const std::size_t Y0 = std::lround(Y);

        // Check that we are inside the domain
        if (X0 < 0 || X0 + 1 >= Width || Y0 < 0 || Y0 + 1 >= Height)
            throw std::runtime_error("Point outside of height map domain.");

        // Compute value by bilinear interpolation
        const double z00 = GridData[Y0 * Width + X0];
        const double z01 = GridData[(Y0 + 1) * Width + X0];
        const double z10 = GridData[Y0 * Width + X0 + 1];
        const double z11 = GridData[(Y0 + 1) * Width + X0 + 1];
        X -= X0;
        Y -= Y0;
        const double z = (1.0 - X) * (1.0 - Y) * z00 + (1.0 - X) * Y * z01 +
                         X * (1.0 - Y) * z10 + X * Y * z11;

        return z;
    }

    // Apply (set) geo reference
    double Apply(const GeoReference& geoReference)
    {
        GridMap = geoReference;
    }

    // Return world (UTM) coordinate of center
    Point2D WorldCoordinateCenter() const
    {
        const Point2D nw = WorldCoordinateNW();
        const Point2D ne = WorldCoordinateNE();
        const Point2D se = WorldCoordinateSE();
        const Point2D sw = WorldCoordinateSW();
        return (nw + ne + se + sw) * 0.25;
    }

    // Return world (UTM) coordinate width
    double WorldCoordinateWidth() const
    {
        const Point2D nw = WorldCoordinateNW();
        const Point2D ne = WorldCoordinateNE();
        return ne.x - nw.x;
    }

    // Return world (UTM) coordinate height
    double WorldCoordinateHeight() const
    {
        const Point2D nw = WorldCoordinateNW();
        const Point2D sw = WorldCoordinateSW();
        return nw.y - sw.y;
    }

    // Return world (UTM) coordinate of NW corner
    Point2D WorldCoordinateNW() const
    {
        return WorldCoordinate(0, 0);
    }

    // Return world (UTM) coordinate of NE corner
    Point2D WorldCoordinateNE() const
    {
        return WorldCoordinate(Width - 1, 0);
    }

    // Return world (UTM) coordinate of SE corner
    Point2D WorldCoordinateSE() const
    {
        return WorldCoordinate(Width - 1, Height - 1);
    }

    // Return world (UTM) coordinate of SW corner
    Point2D WorldCoordinateSW() const
    {
        return WorldCoordinate(0, Height - 1);
    }

    // Return world (UTM) coordinate of pixel (X, Y) = (column, row)
    Point2D WorldCoordinate(size_t X, size_t Y) const
    {
        const double AX = GridMap.A * X;
        const double BY = GridMap.B * Y;
        const double DX = GridMap.D * X;
        const double EY = GridMap.E * Y;
        const double C = GridMap.C;
        const double F = GridMap.F;
        Point2D p(AX + BY + C, DX + EY + F);
        return p;
    }

};

std::ostream& operator<<(std::ostream& stream, const HeightMap& heightMap)
{
    stream << "C = " << heightMap.WorldCoordinateCenter()
           << " W = " << heightMap.WorldCoordinateWidth()
           << " H = " << heightMap.WorldCoordinateHeight()
           << " NW = " << heightMap.WorldCoordinateNW()
           << " NE = " << heightMap.WorldCoordinateNE()
           << " SE = " << heightMap.WorldCoordinateSE()
           << " SW = " << heightMap.WorldCoordinateSW();
}

}

#endif
