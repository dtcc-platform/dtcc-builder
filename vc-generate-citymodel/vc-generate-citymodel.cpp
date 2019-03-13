// VirtualCity@Chalmers: vc-generate-citymodel
// Anders Logg 2019

#include <iostream>
#include <iomanip>
#include <random>

#include "CityModel.h"
#include "OSM.h"
#include "JSON.h"

// FIXME: Testing
#include "CoordinateSystem.h"

using namespace VirtualCity;

void Help()
{
    std::cerr << "Usage: vc-generate-citymodel OpenStreetMap.osm Parameters.json"
              << std::endl;
}

CityModel GenerateCityModel(std::string fileNameOpenStreetMap)
{
    // Create empty city model
    CityModel cityModel;

    // Read OSM data from file
    OSM::Read(cityModel, fileNameOpenStreetMap);

    // FIXME: Temporary until we have the building heights
    for (size_t i = 0; i < cityModel.Buildings.size(); i++)
        cityModel.Buildings[i].Height = 30.0;

    return cityModel;
}

int main(int argc, char* argv[])
{
    // Check command-line arguments
    if (argc != 3)
    {
        Help();
        return 1;
    }

    // Get filenames
    const std::string fileNameOpenStreetMap(argv[1]);
    const std::string fileNameParameters(argv[2]);

    // Read parameters from file
    Parameters parameters;
    JSON::Read(parameters, fileNameParameters);

    // Report used parameters
    // FIXME:

    // Generate city model
    CityModel cityModel = GenerateCityModel(fileNameOpenStreetMap);
    std::cout << cityModel << std::endl;

    // Write to file
    JSON::Write(cityModel, "CityModel.json");

    // FIXME: Testing 4326
    Point2D p(0, 0);
    Point2D q = CoordinateSystem::Transform(p, "epsg:4326", "epsg:3007");
    Point2D r = CoordinateSystem::Transform(p, "epsg:3006", "epsg:3007");
    std::cout << p << std::endl;
    std::cout << q << std::endl;
    std::cout << r << std::endl;

    return 0;
}
