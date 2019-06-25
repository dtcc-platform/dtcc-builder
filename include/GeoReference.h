// Georeference following the Esri World File format.
// https://en.wikipedia.org/wiki/World_file
// Anders Logg 2019

#ifndef VC_GEO_REFERENCE_H
#define VC_GEO_REFERENCE_H

#include <pair>
#include <string>

namespace VirtualCity
{

class GeoReference
{
public:
  double A; // x-component of pixel width
  double D; // y-component of pixel width
  double B; // x-component of pixel height
  double E; // y-component of pixel height
  double C; // x-coordinate of center of upper left pixel
  double F; // y-coordinate of center of upper left pixel

  // Create empty (zero) geo reference
  GeoReference() : A(0), D(0), B(0), E(0), C(0), F(0) {}

  // Map pixel to world (UTM) coordinates
  Point2D Pixel2Coordinate(int X, int Y) const
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

  // Map world (UTM) coordinates to pixels
  std::pair<int, int> Coordinate2Pixel(const Point2D &p) const
  {
    const double Ex = GridMap.E * p.x;
    const double By = GridMap.B * p.y;
    const double BF = GridMap.B * GridMap.F;
    const double EC = GridMap.E * GridMap.C;
    const double Dx = GridMap.D * p.x;
    const double Ay = GridMap.A * p.y;
    const double DC = GridMap.D * GridMap.C;
    const double AF = GridMap.A * GridMap.F;
    const double AE = GridMap.A * GridMap.E;
    const double DB = GridMap.D * GridMap.B;
    const double det = AE - DB;
    double X = (Ex - By + BF - EC) / det;  // column
    double Y = (-Dx + Ay + DC - AF) / det; // row

    std::pair<int, int> pixel;
    pixel.first = std::lround(X);
    pixel.second = std::lround(Y);

    return pixel;
  }
};

} // namespace VirtualCity

#endif
