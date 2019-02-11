// Georeference following the Esri World File format.
// https://en.wikipedia.org/wiki/World_file
// Anders Logg 2019

#ifndef VC_GEO_REFERENCE_H
#define VC_GEO_REFERENCE_H

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

};

}

#endif
