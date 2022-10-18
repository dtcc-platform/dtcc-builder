// STL I/O
// Anders Logg 2018
// Licensed under the MIT License

#ifndef DTCC_STL_H
#define DTCC_STL_H

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "Logging.h"

namespace DTCC_BUILDER
{

class STL
{
public:
  // Read (binary) STL file and return triangles as a flattened
  // list of vertex coordinates (3 x 3 = 9 numbers per triangle).
  static std::vector<float> Read(const std::string& filename)
  {
    info("STL file:   " + filename);

    // Open file
    std::ifstream f(filename.c_str(), std::ios::in | std::ios::binary);

    // Check file
    if (!f)
      throw std::runtime_error("Unable to open file: " + filename);

    // Parse header and number of triangles
    std::string InfoHeader = ParseHeader(f);
    unsigned long NumTriangles = ParseUnsignedLong(f);
    info("STL header: " + InfoHeader);
    info("STL size:  " + str(NumTriangles) + " triangles");

    // Initialize mesh data
    std::vector<float> Triangles(9 * NumTriangles);

    // Parse data
    for (unsigned long i = 0; i < NumTriangles; i++)
    {
      // Parse normal (not used)
      float nx = ParseFloat(f);
      float ny = ParseFloat(f);
      float nz = ParseFloat(f);

      // Parse triangle vertices
      for (int j = 0; j < 9; j++)
        Triangles[9 * i + j] = ParseFloat(f);

      // Parse dummy 2 bytes
      char buf[2];
      f.read(buf, 2);
    }

    return Triangles;
  }

private:
  static std::string ParseHeader(std::ifstream &f)
  {
    char InfoHeader[80] = "";
    f.read(InfoHeader, 80);
    for (int i = 0; i < 80; i++)
      if (InfoHeader[i] == '\n')
        InfoHeader[i] = '\0'; // cleanup
    return std::string(InfoHeader);
  }

  static unsigned long ParseUnsignedLong(std::ifstream &f)
  {
    char buf[4];
    f.read(buf, 4);
    return *((unsigned long *)buf);
  }

  static float ParseFloat(std::ifstream &f)
  {
    char buf[sizeof(float)];
    f.read(buf, 4);
    return *((float *)buf);
  }
};

} // namespace DTCC_BUILDER

#endif
