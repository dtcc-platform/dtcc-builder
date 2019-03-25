// VirtualCity@Chalmers: vc-generate-citymodel
// Anders Logg 2019

#include <iostream>
#include <iomanip>
#include <random>

#include "CommandLine.h"
#include "OSM.h"
#include "SHP.h"
#include "JSON.h"
#include "CityModel.h"
#include "CityModelGenerator.h"

using namespace VirtualCity;

void Help()
{
    std::cerr << "Usage: vc-generate-citymodel PropertyMap.[osm/shp]"
              << " HeightMap.json" << std::endl;
}

int main(int argc, char* argv[])
{
    // Check command-line arguments
    if (argc != 4)
    {
        Help();
        return 1;
    }

    // Get filenames
    const std::string fileNamePropertyMap(argv[1]);
    const std::string fileNameHeightMap(argv[2]);

    // Read polygons
    std::vector<Polygon> polygons;
    if (CommandLine::EndsWith(fileNamePropertyMap, ".osm"))
        OSM::Read(polygons, fileNamePropertyMap);
    else if (CommandLine::EndsWith(fileNamePropertyMap, ".shp"))
        SHP::Read(polygons, fileNamePropertyMap);

    // Read height map
    HeightMap heightMap;
    JSON::Read(heightMap, fileNameHeightMap);

    // Generate city model
    CityModel cityModel;
    CityModelGenerator::GenerateCityModel(cityModel, polygons, heightMap);

    // Pretty-print
    std::cout << cityModel << std::endl;

    // Write to file
    JSON::Write(cityModel, "CityModel.json");

    return 0;
}
