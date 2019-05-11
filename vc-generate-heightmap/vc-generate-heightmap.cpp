// VirtualCity@Chalmers: vc-generate-heightmap
// Anders Logg 2019

#include <iostream>
#include <string>
#include <vector>

#include "CommandLine.h"
#include "HeightMap.h"
#include "HeightMapGenerator.h"
#include "LAS.h"
#include "JSON.h"

using namespace VirtualCity;

void help()
{
    std::cerr << "Usage: vc-generate-heightmap "
              << "PointCloud0.las PointCloud1.las ... Parameters.json"
              << std::endl;
}

int main(int argc, char* argv[])
{
    // Check command-line arguments
    if (argc < 3)
    {
        help();
        return 1;
    }

    // Get filenames
    std::vector<std::string> fileNamesLAS;
    for (size_t i = 1; i < argc - 1; i++)
        fileNamesLAS.push_back(std::string(argv[i]));
    const std::string fileNameParameters(argv[argc - 1]);

    // Read parameters from file
    Parameters parameters;
    JSON::Read(parameters, fileNameParameters);
    std::cout << parameters << std::endl;

    // Read point cloud from LAS files
    PointCloud pointCloud;
    for (const std::string fileNameLAS : fileNamesLAS)
    {
        LAS::Read(pointCloud, fileNameLAS);
        std::cout << pointCloud << std::endl;
    }

    // Generate height map
    HeightMap heightMap;
    HeightMapGenerator::GenerateHeightMap(heightMap, pointCloud,
                                          parameters.X0, parameters.Y0,
                                          parameters.XMin, parameters.YMin,
                                          parameters.XMax, parameters.YMax,
                                          parameters.HeightMapResolution);
    std::cout << heightMap << std::endl;

    // Write height map to JSON file
    JSON::Write(heightMap, "HeightMap.json");

    return 0;
}
