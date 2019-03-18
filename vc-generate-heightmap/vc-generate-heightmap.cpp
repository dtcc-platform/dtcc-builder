// VirtualCity@Chalmers: vc-generate-heightmap
// Anders Logg 2019

#include <iostream>
#include <string>

#include "CommandLine.h"
#include "HeightMap.h"
#include "HeightMapGenerator.h"
#include "LAS.h"
#include "JSON.h"

using namespace VirtualCity;

void help()
{
    std::cerr << "Usage: vc-generate-heightmap "
              << "PointCloud.las Parameters.json"
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
    const std::string fileNameLAS(argv[1]);
    const std::string fileNameParameters(argv[2]);

    // Read parameters from file
    Parameters parameters;
    JSON::Read(parameters, fileNameParameters);

    // Report used parameters
    std::cout << "vc-generate-heightmap: XMin = "
              << parameters.XMin << std::endl;
    std::cout << "vc-generate-heightmap: YMin = "
              << parameters.YMin << std::endl;
    std::cout << "vc-generate-heightmap: XMax = "
              << parameters.XMax << std::endl;
    std::cout << "vc-generate-heightmap: YMin = "
              << parameters.YMax << std::endl;
    std::cout << "vc-generate-heightmap: HeightMapResolution = "
              << parameters.HeightMapResolution << std::endl;

    // Read point cloud from LAS file
    PointCloud pointCloud;
    LAS::Read(pointCloud, fileNameLAS);

    // Generate height map
    HeightMap heightMap(parameters.XMin, parameters.YMin,
                        parameters.XMax, parameters.YMax,
                        parameters.HeightMapResolution);
    HeightMapGenerator::GenerateHeightMap(heightMap, pointCloud);

    // Write height map to JSON file
    JSON::Write(heightMap, "HeightMap.json");

    return 0;
}
