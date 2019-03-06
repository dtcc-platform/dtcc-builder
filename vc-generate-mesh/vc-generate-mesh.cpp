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
#include "FEniCS.h"

#include "CSV.h"

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

    // Set filename for output
    size_t idx = fileNameCityModel.rfind('.');
    const std::string fileNamePrefix = fileNameCityModel.substr(0, idx);
    std::cout << "Prefix: " << fileNamePrefix << std::endl;

    // Read city model from file
    CityModel cityModel;
    JSON::Read(cityModel, fileNameCityModel);

    // Read height map from file
    HeightMap heightMap;
    JSON::Read(heightMap, fileNameHeightMap);
    std::cout << heightMap << std::endl;

    // Read parameters from file
    Parameters parameters;
    JSON::Read(parameters, fileNameParameters);

    // Report used parameters
    const Point2D C(parameters.DomainCenterX, parameters.DomainCenterY);
    const double R = parameters.DomainRadius;
    const double H = parameters.DomainHeight;
    const double h = parameters.MeshSize;
    std::cout << "vc-generate-mesh: DomainCenter = " << C << std::endl;
    std::cout << "vc-generate-mesh: DomainRadius = " << R << std::endl;
    std::cout << "vc-generate-mesh: DomainHeight = " << H << std::endl;
    std::cout << "vc-generate-mesh: MeshSize = " << h << std::endl;

    // Generate 2D mesh
    Mesh2D mesh2D = MeshGenerator::GenerateMesh2D(cityModel, C, R, h);

    // Generate mesh (excluding height map)
    Mesh3D mesh3D = MeshGenerator::GenerateMesh3D(mesh2D, cityModel, H, h);

    // Convert to FEniCS meshes
    dolfin::Mesh _mesh2D, _mesh3D;
    FEniCS::ConvertMesh(mesh2D, _mesh2D);
    FEniCS::ConvertMesh(mesh3D, _mesh3D);

    // Apply mesh smoothing to account for height map
    MeshSmoother::SmoothMesh(_mesh3D,
                             heightMap,
                             cityModel,
                             mesh3D.DomainMarkers,
                             H, h);

    // Generate height map function (used only for testing/visualization)
    auto z = MeshSmoother::GenerateHeightMapFunction(_mesh2D, heightMap);

    // Generate mesh boundary (used only for testing/visualization)
    dolfin::BoundaryMesh _boundary3D(_mesh3D, "exterior");

    // Write to files
    dolfin::File f0(fileNamePrefix + "Mesh.xml");
    dolfin::File f1(fileNamePrefix + "Mesh.pvd");
    dolfin::File f2(fileNamePrefix + "Boundary.pvd");
    dolfin::File f3(fileNamePrefix + "HeightMap.pvd");
    f0 << _mesh3D;
    f1 << _mesh3D;
    f2 << _boundary3D;
    f3 << *z;

    return 0;
}
