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

    // Get filename
    string filename(argv[1]);

    // FIXME: Get filename from command-line arguments
    std::string fileName = "HeightMap.json";

    // Read height map from file
    HeightMap heightmap;
    JSON::Read(heightmap, fileName);

    // Generate mesh (excluding height map)
    Mesh3D m = MeshGenerator::GenerateMesh3D();

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
