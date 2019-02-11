// PNG I/O
// Anders Logg 2019

#ifndef VC_PNG_H
#define VC_PNG_H

#include <stdexcept>
#include <iostream>
#include <stdio.h>
#include <Magick++.h>

#include "HeightMap.h"

namespace VirtualCity
{

class PNG
{
public:

    // Read height map from PNG file
    static void Read(HeightMap& heightMap, std::string fileName)
    {
        std::cout << "PNG: " << "Reading height map from file "
                  << fileName << std::endl;

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
        std::cout << std::fixed
                  << "PNG: " << numCols << " x " << numRows
                  << " (" << numPixels << " pixels)" << std::endl;

        // Write pixel data to array
        std::vector<int16_t> pixels(numPixels);
        image.write(0, 0, numCols, numRows, "R",
                    Magick::ShortPixel, pixels.data());

        // Convert to floating point height map
        heightMap.Width = numCols;
        heightMap.Height = numRows;
        heightMap.GridData.resize(numPixels);
        for (size_t i = 0; i < numPixels; i++)
        {
            double z = pixels[i];
            z -= 1024; // 1024 is sea level
            z *= 0.01; // pixel data is cm
            heightMap.GridData[i] = z;
        }
    };

};

}

#endif
