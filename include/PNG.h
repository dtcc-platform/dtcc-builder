// PNG I/O
// Anders Logg 2019
// Licensed under the MIT License

#ifndef DTCC_PNG_H
#define DTCC_PNG_H

#include <Magick++.h>
#include <iostream>
#include <stdexcept>
#include <stdio.h>

#include "HeightMap.h"
#include "Logging.h"

namespace DTCC
{

class PNG
{
public:
  // Read height map from PNG file
  static void Read(HeightMap &heightMap, std::string fileName, size_t stride)
  {
    Info("PNG: Reading height map from file " + fileName);

    // Read image
    Magick::Image image;
    image.read(fileName);

    // Check type (should be grayscale)
    if (image.type() != Magick::GrayscaleType)
      throw std::runtime_error("Illegal image type; expecting grayscale.");

    // Check depth (should be 16)
    if (image.depth() != 16)
      throw std::runtime_error("Illegal image depth; expecting 16.");

    // Get image dimensions
    const size_t numRows = image.rows();
    const size_t numCols = image.columns();
    const size_t numPixels = numRows * numCols;
    Info("PNG: " + str(numCols) + " x " + str(numRows) + " (" + str(numPixels) + " pixels)");

    // Write pixel data to array
    std::vector<double> pixels(numPixels);
    image.write(0, 0, numCols, numRows, "I", Magick::DoublePixel,
                pixels.data());

    // Convert to floating point height map using stride
    heightMap.Width = (numCols - 1) / stride + 1;
    heightMap.Height = (numRows - 1) / stride + 1;
    const size_t numScaledPixels = heightMap.Width * heightMap.Height;
    heightMap.GridData.resize(numScaledPixels);
    size_t k = 0;
    for (size_t i = 0; i < numRows; i += stride)
    {
      for (size_t j = 0; j < numCols; j += stride)
      {
        double z = pixels[i * numCols + j];
        z *= 65535; // 16 bit depth
        z -= 1024;  // 1024 is sea level
        z *= 0.01;  // pixel data is cm
        Progress("z = " + z);
        heightMap.GridData[k++] = z;
      }
    }
  }
};

} // namespace DTCC

#endif
