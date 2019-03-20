// VirtualCity@Chalmers: vc-generate-citymodel
// Anders Logg 2019

#include <iostream>
#include <iomanip>
#include <random>

#include "CommandLine.h"
#include "CityModel.h"
#include "OSM.h"
#include "SHP.h"
#include "JSON.h"

// FIXME: Testing
#include "CoordinateSystem.h"

using namespace VirtualCity;

void Help()
{
    std::cerr << "Usage: vc-generate-citymodel PropertyMap.[osm/shp]"
              << " Parameters.json" << std::endl;
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
    const std::string fileNameFootprints(argv[1]);
    const std::string fileNameParameters(argv[2]);

    // Read parameters from file
    Parameters parameters;
    JSON::Read(parameters, fileNameParameters);

    // Report used parameters
    // FIXME: Not implemented

    // Extract footprints from property map
    std::vector<Polygon> polygons;
    if (CommandLine::EndsWith(fileNameFootprints, ".osm"))
        OSM::Read(polygons, fileNameFootprints);
    else if (CommandLine::EndsWith(fileNameFootprints, ".shp"))
        SHP::Read(polygons, fileNameFootprints);

    CityModel cityModel;
    std::cout << cityModel << std::endl;

    // Extract heights
    // FIXME: Not implemented

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
