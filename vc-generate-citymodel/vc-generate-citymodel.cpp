// VirtualCity@Chalmers: vc-generate-citymodel
// Anders Logg 2019

#include <iostream>

#include "CityModel.h"
#include "JSON.h"

using namespace VirtualCity;

void help()
{
    std::cerr << "Usage: vc-generate-citymodel"
              << std::endl;
}

int main(int argc, char* argv[])
{
    // Check command-line arguments
    if (argc != 1)
    {
        help();
        return 1;
    }

    // FIXME: This is just a test script for now

    Building b0;
    b0.Footprint.push_back(Point2D(-2, -2));
    b0.Footprint.push_back(Point2D(-1, -2));
    b0.Footprint.push_back(Point2D(-1, -1));
    b0.Footprint.push_back(Point2D(-2, -1));

    Building b1;
    b1.Footprint.push_back(Point2D(2, 3));
    b1.Footprint.push_back(Point2D(3, 3));
    b1.Footprint.push_back(Point2D(3, 4));
    b1.Footprint.push_back(Point2D(2, 4));

    CityModel cityModel;
    cityModel.Buildings.push_back(b0);
    cityModel.Buildings.push_back(b1);

    JSON::Write(cityModel, "CityModel.json");

    CityModel test;
    JSON::Read(test, "CityModel.json");
    JSON::Write(test, "test.json");

    return 0;
}
