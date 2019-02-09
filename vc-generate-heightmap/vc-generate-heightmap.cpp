// VirtualCity@Chalmers: vc-generate-heightmap
// Anders Logg 2019

#include <iostream>

#include "JSON.h"
#include "PNG.h"

using namespace VirtualCity;

/*
char* getopt(std::string option, int argc, char* argv[])
{
    auto end = argc + argv;
    auto it = find(argv, end, option);
    if (it != end && ++it != end)
        return *it;
    return 0;
}
*/

void help()
{
    std::cerr << "Usage: vc-generate-heightmap [options] HeightMap.png" << std::endl;
}

int main(int argc, char* argv[])
{
    /*

    // Check command-line arguments
    if (argc != 2)
    {
        help();
        return 1;
    }

    // Get filename
    string filename(argv[1]);

    */

    // FIXME: Test data

    // Read height map from PNG file
    HeightMap heightMap;
    PNG::Read(heightMap, "HeightMapLindholmen.png");

    // Write height map to JSON file
    JSON::Write(heightMap, "HeightMapLindholmen.json");

    return 0;
}
