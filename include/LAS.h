// LAS I/O
// Anders Logg 2019

#ifndef VC_SHP_H
#define VC_SHP_H

#include <iostream>
#include <fstream>
#include <string>
#include <liblas/liblas.hpp>

#include "SurfaceModel.h"
#include "Point.h"

namespace VirtualCity
{

class LAS
{
public:

    // Read height map from LAS file
    static void Read(SurfaceModel& surfaceModel,
                     std::string fileName)
    {
        std::cout << "LAS: " << "Reading surface model from file "
                  << fileName << std::endl;

        // Open file
        std::ifstream f;
        f.open(fileName, std::ios::in | std::ios::binary);

        // Create reader
        liblas::ReaderFactory factory;
        liblas::Reader reader = factory.CreateWithStream(f);

        // Read header
        liblas::Header const& header = reader.GetHeader();
        const bool isCompressed = header.Compressed();
        const std::string signature = header.GetFileSignature();
        const size_t numPoints = header.GetPointRecordsCount();
        if (isCompressed)
            std::cout << "LAS: Compressed" << std::endl;
        else
            std::cout << "LAS: Uncompressed" << std::endl;
        std::cout << "LAS: " << signature << std::endl;
        std::cout << "LAS: " << numPoints << " points" << std::endl;

        // Iterate over points
        while (reader.ReadNextPoint())
        {
            // Get point
            liblas::Point const& _p = reader.GetPoint();
            const Point3D p(_p.GetX(), _p.GetY(), _p.GetZ());

            std::cout << p << std::endl;

        }

    }

};

}

#endif
