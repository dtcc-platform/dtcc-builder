// VirtualCity@Chalmers: vc-generate-mesh
// Anders Logg 2018

#include <iostream>
#include <algorithm>

// FIXME: Should not be used here
#include <dolfin.h>

#include "CommandLine.h"
#include "Parameters.h"
#include "MeshGenerator.h"
#include "MeshSmoother.h"
#include "HeightMap.h"
#include "JSON.h"

using namespace std;
using namespace VirtualCity;

void help()
{
    cerr << "Usage: vc-generate-mesh CityModel.json HeightMap.json Parameters.json" << endl;
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
    const std::string fileNameCityModel(argv[1]);
    const std::string fileNameHeightMap(argv[2]);
    const std::string fileNameParameters(argv[3]);

    // Read city model from file
    CityModel cityModel;
    JSON::Read(cityModel, fileNameCityModel);

    // Read height map from file
    HeightMap heightMap;
    JSON::Read(heightMap, fileNameHeightMap);

    // Read parameters from file
    Parameters parameters;
    JSON::Read(parameters, fileNameParameters);

    // Report used parameters
    std::cout << "vc-generate-mesh: DomainRadius = "
              << parameters.DomainRadius << std::endl;
    std::cout << "vc-generate-mesh: MeshSize = "
              << parameters.MeshSize << std::endl;

    // Generate mesh (excluding height map)
    Mesh3D m = MeshGenerator::GenerateMesh3D(cityModel);

    // Apply mesh smoothing to account for height map
    MeshSmoother::SmoothMesh();

    // FIXME: Write test output
    dolfin::Mesh _m("Mesh3D.xml");
    dolfin::BoundaryMesh _b(_m, "exterior");
    dolfin::File _f("Mesh3D.pvd");
    dolfin::File _g("Mesh3DBoundary.pvd");
    _f << _m;
    _g << _b;

    return 0;
}
