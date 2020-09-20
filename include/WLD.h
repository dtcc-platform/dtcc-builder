// WLD I/O
// Anders Logg 2019
// Licensed under the MIT License

#ifndef DTCC_WLD_H
#define DTCC_WLD_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "GeoReference.h"
#include "Parameters.h"

namespace DTCC
{

class WLD
{
public:
  static const int PRECISION = 16;

  // Read geo reference from WLD file
  static void
  Read(GeoReference &geoReference, std::string fileName, size_t stride = 1)
  {
    std::cout << "WLD: "
              << "Reading geo reference from file " << fileName << std::endl;

    std::ifstream f(fileName);
    std::string A, D, B, E, C, F;

    f >> A;
    f >> D;
    f >> B;
    f >> E;
    f >> C;
    f >> F;

    geoReference.A = std::stod(A) * double(stride);
    geoReference.D = std::stod(D) * double(stride);
    geoReference.B = std::stod(B) * double(stride);
    geoReference.E = std::stod(E) * double(stride);
    geoReference.C = std::stod(C);
    geoReference.F = std::stod(F);
  };

  // Write geo reference to WLD file
  static void Write(GeoReference &geoReference, std::string fileName)
  {
    std::cout << "WLD: "
              << "Writing geo reference to file " << fileName << std::endl;

    std::ofstream f(fileName);
    f << std::setprecision(Parameters::Precision);

    f << geoReference.A << std::endl;
    f << geoReference.D << std::endl;
    f << geoReference.B << std::endl;
    f << geoReference.E << std::endl;
    f << geoReference.C << std::endl;
    f << geoReference.F << std::endl;
  };
};

} // namespace DTCC

#endif
