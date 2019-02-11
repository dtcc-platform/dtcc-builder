// Representation of a 2.5D city model.
// Copyright (C) 2019 Anders Logg.

#ifndef VC_CITY_MODEL_H
#define VC_CITY_MODEL_H

#include <vector>

#include "Building.h"

namespace VirtualCity
{

class CityModel
{
public:

    // List of buildings
    std::vector<Building> Buildings;

};

}

#endif
