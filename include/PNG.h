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

        // Get image dimensions
        const size_t numRows = image.rows();
        const size_t numCols = image.columns();

        std::cout << "rows = " << numRows << std::endl;
        std::cout << "cols = " << numCols << std::endl;

        for (size_t i = 0; i < numRows; i++)
        {
            for (size_t j = 0; j < numCols; j++)
            {
                Magick::Color color = image.pixelColor(j, i);
                std::cout << i << " " << j << " " << color.redQuantum() << std::endl;
            }
        }

    };

};

}

#endif
