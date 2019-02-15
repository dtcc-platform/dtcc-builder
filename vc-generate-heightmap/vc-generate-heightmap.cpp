// VirtualCity@Chalmers: vc-generate-heightmap
// Anders Logg 2019

#include <iostream>
#include <string>

#include "CommandLine.h"
#include "HeightMap.h"
#include "GeoReference.h"
#include "PNG.h"
#include "WLD.h"
#include "JSON.h"

using namespace VirtualCity;

void help()
{
    std::cerr << "Usage: vc-generate-heightmap "
              << "HeightMap.png HeightMap.wld HeightMap.json Parameters.json"
              << std::endl;
}

int main(int argc, char* argv[])
{
    // Check command-line arguments
    if (argc < 5)
    {
        help();
        return 1;
    }

    // Get filenames
    const std::string fileNamePNG(argv[1]);
    const std::string fileNameWLD(argv[2]);
    const std::string fileNameJSON(argv[3]);
    const std::string fileNameParameters(argv[4]);

    // Read parameters from file
    Parameters parameters;
    JSON::Read(parameters, fileNameParameters);

    // Report used parameters
    std::cout << "vc-generate-mesh: HeightMapStride = "
              << parameters.HeightMapStride << std::endl;

    // Read height map from PNG file
    HeightMap heightMap;
    PNG::Read(heightMap, fileNamePNG, parameters.HeightMapStride);

    // Read geo reference from WLD file
    GeoReference geoReference;
    WLD::Read(geoReference, fileNameWLD, parameters.HeightMapStride);

    // Apply geo reference to height map
    heightMap.Apply(geoReference);
    std::cout << heightMap << std::endl;

    // Write height map to JSON file
    JSON::Write(heightMap, fileNameJSON);

    return 0;
}
