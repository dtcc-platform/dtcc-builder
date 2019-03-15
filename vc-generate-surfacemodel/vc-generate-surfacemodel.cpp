// VirtualCity@Chalmers: vc-generate-heightmap
// Anders Logg 2019

#include <iostream>
#include <string>

#include "CommandLine.h"
#include "HeightMap.h"
#include "LAS.h"
#include "JSON.h"

using namespace VirtualCity;

void help()
{
    std::cerr << "Usage: vc-generate-surfacemodel "
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
    std::cout << "vc-generate-mesh: HeightMapStride = "
              << parameters.HeightMapStride << std::endl;

    // Read point cloud from LAS file
    PointCloud pointCloud;
    LAS::Read(pointCloud, fileNameLAS);

    // // Read geo reference from WLD file
    // GeoReference geoReference;
    // WLD::Read(geoReference, fileNameWLD, parameters.HeightMapStride);

    // // Apply geo reference to height map
    // heightMap.Apply(geoReference);
    // std::cout << heightMap << std::endl;

    // Write height map to JSON file
    //JSON::Write(surfaceModel, "SurfaceModel.json");

    return 0;
}
