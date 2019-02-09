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

/* TODO

[x] remove using namespace std
[x] plotting utilty
[x] link to triangle
[x] parse command-line
[x] generate convex hull
[x] generate 2D mesh
[x] generate 3D mesh
[x] save to FEniCS format
[ ] timings
[x] handle buildings when generating layers
[ ] handle height map
[ ] parsing of building and height-map data
[ ] mesh smoothing
[ ] don't append .csv etc in I/O
[ ] unified interface Read/Write for I/O
[ ] use C# coding style

*/

void help()
{
    cerr << "Usage: vc-generate-mesh [options] mesh.stl" << endl;
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

    // FIXME: Get filename from command-line arguments
    std::string fileName = "HeightMap.json";

    // Read height map from file
    HeightMap heightmap;
    //JSON::Read(heightmap, fileName);

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
