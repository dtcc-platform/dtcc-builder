// PNG I/O
// Anders Logg 2019

#ifndef PNG_H
#define PNG_H

#include <iostream>
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

        std::cout << "Nothing to see here yet..." << std::endl;
    };

};

}

#endif
