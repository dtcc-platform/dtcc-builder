// VirtualCity@Chalmers: vc-generate-heightmap
// Anders Logg 2019

#include <iostream>
#include <string>

#include "HeightMap.h"
#include "GeoReference.h"
#include "PNG.h"
#include "WLD.h"
#include "JSON.h"

using namespace VirtualCity;

void help()
{
    std::cerr << "Usage: vc-generate-heightmap "
              << "HeightMap.png HeightMap.wld HeightMap.json"
              << std::endl;
}

int main(int argc, char* argv[])
{
    // Check command-line arguments
    if (argc != 4)
    {
        help();
        return 1;
    }

    // Get filenames
    const std::string fileNamePNG(argv[1]);
    const std::string fileNameWLD(argv[2]);
    const std::string fileNameJSON(argv[3]);

    // Read height map from PNG file
    HeightMap heightMap;
    PNG::Read(heightMap, fileNamePNG);

    // Read geo reference from WLD file
    GeoReference geoReference;
    WLD::Read(geoReference, fileNameWLD);

    // Apply geo reference to height map
    heightMap.Apply(geoReference);

    // Write height map to JSON file
    JSON::Write(heightMap, fileNameJSON);

    return 0;
}
