// VirtualCity@Chalmers: vc-generate-citymodel
// Anders Logg 2019

#include <iostream>
#include <iomanip>
#include <random>

#include "CityModel.h"
#include "OSM.h"
#include "JSON.h"

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

    // Write to file
    JSON::Write(cityModel, "CityModel.json");//

    return 0;
}
