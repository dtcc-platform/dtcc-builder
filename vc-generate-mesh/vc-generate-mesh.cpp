// VirtualCity@Chalmers: vc-generate-mesh
// Anders Logg 2018

#include <iostream>
#include <algorithm>

// FIXME: Should not be used here
#include <dolfin.h>

#include "CommandLine.h"
#include "MeshGenerator.h"
#include "MeshSmoother.h"
#include "HeightMap.h"
#include "JSON.h"

using namespace std;
using namespace VirtualCity;

void help()
{
    cerr << "Usage: vc-generate-mesh CityModel.json HeightMap.json" << endl;
}

int main(int argc, char* argv[])
{
    // Check command-line arguments
    if (argc != 3)
    {
        help();
        return 1;
    }

    // Get filenames
    const std::string fileNameCityModel(argv[1]);
    const std::string fileNameHeightMap(argv[2]);

    // Read city model from file
    CityModel cityModel;
    JSON::Read(cityModel, fileNameCityModel);

    // Read height map from file
    HeightMap heightMap;
    JSON::Read(heightMap, fileNameHeightMap);

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
