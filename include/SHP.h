// SHP I/O
// Anders Logg 2019

#ifndef VC_SHP_H
#define VC_SHP_H

#include <iostream>
#include <shapefil.h>

namespace VirtualCity
{

class SHP
{
public:

    // Read city model from SHP file. Note that the corresponding
    // .shx and .dbf files must also be present in the same directory.
    static void Read(CityModel& cityModel, std::string fileName)
    {
        SHPOpen(fileName.c_str(), "r");

    }

};

}

#endif
