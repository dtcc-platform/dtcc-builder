// VirtualCity@Chalmers: vc-generate-citymodel
// Anders Logg 2019

#include <iostream>
#include <iomanip>
#include <random>

#include "CityModel.h"
#include "Geometry.h"
#include "JSON.h"

using namespace VirtualCity;

void Help()
{
    std::cerr << "Usage: vc-generate-citymodel Parameters.json"
              << std::endl;
}

CityModel GenerateCityModel()
{
    // Create empty city model
    CityModel cityModel;


    return cityModel;
}

int main(int argc, char* argv[])
{
    // Check command-line arguments
    if (argc != 2)
    {
        Help();
        return 1;
    }

    // Get filenames
    const std::string fileNameParameters(argv[1]);

    // Read parameters from file
    Parameters parameters;
    JSON::Read(parameters, fileNameParameters);

    // Report used parameters
    // FIXME:

    // Generate city model
    CityModel cityModel = GenerateCityModel();

    // Write to file
    JSON::Write(cityModel, "CityModel.json");//

    return 0;
}
