// Georeference following the Esri World File format.
// https://en.wikipedia.org/wiki/World_file
// Anders Logg 2019

#ifndef DTCC_GEO_REFERENCE_H
#define DTCC_GEO_REFERENCE_H

#include <utility>
	
// #include <pair>
#include <string>

namespace DTCC
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
  Vector2D Pixel2Coordinate(int X, int Y) const
  {
    const double AX = A * X;
    const double BY = B * Y;
    const double DX = D * X;
    const double EY = E * Y;

    Vector2D p(AX + BY + C, DX + EY + F);
    return p;
  }

  // Map world (UTM) coordinates to pixels
  std::pair<int, int> Coordinate2Pixel(const Vector2D &p) const
  {
    const double Ex = E * p.x;
    const double By = B * p.y;
    const double BF = B * F;
    const double EC = E * C;
    const double Dx = D * p.x;
    const double Ay = A * p.y;
    const double DC = D * C;
    const double AF = A * F;
    const double AE = A * E;
    const double DB = D * B;
    const double det = AE - DB;
    double X = (Ex - By + BF - EC) / det;  // column
    double Y = (-Dx + Ay + DC - AF) / det; // row

    std::pair<int, int> pixel;
    pixel.first = std::lround(X);
    pixel.second = std::lround(Y);

    return pixel;
  }
};

} // namespace DTCC

#endif
