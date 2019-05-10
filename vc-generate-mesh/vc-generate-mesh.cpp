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

    // Read city model from file
    CityModel cityModel;
    JSON::Read(cityModel, fileNameCityModel);
    std::cout << cityModel << std::endl;

    // Read height map from file
    HeightMap heightMap;
    JSON::Read(heightMap, fileNameHeightMap);
    std::cout << heightMap << std::endl;

    // Read parameters from file
    Parameters parameters;
    JSON::Read(parameters, fileNameParameters);

    // Report used parameters
    std::cout << "vc-generate-mesh: MeshResolution = "
              << parameters.MeshResolution << std::endl;
    std::cout << "vc-generate-mesh: DomainHeight = "
              << parameters.DomainHeight << std::endl;

    // Generate 2D mesh
    Mesh2D mesh2D = MeshGenerator::GenerateMesh2D(cityModel,
                    heightMap.XMin,
                    heightMap.YMin,
                    heightMap.XMax,
                    heightMap.YMax,
                    parameters.MeshResolution);
    std::cout << mesh2D << std::endl;

    // Generate mesh (excluding height map)
    // Mesh3D mesh3D = MeshGenerator::GenerateMesh3D(mesh2D,
    //                 cityModel,
    //                 parameters.DomainHeight,
    //                 parameters.MeshResolution);
    // std::cout << mesh3D << std::endl;

    // Convert to FEniCS meshes
    dolfin::Mesh _mesh2D, _mesh3D;
    FEniCS::ConvertMesh(mesh2D, _mesh2D);
    //FEniCS::ConvertMesh(mesh3D, _mesh3D);

    // FIXME: Testing
    // Apply mesh smoothing to account for height map
    // MeshSmoother::SmoothMesh(_mesh3D,
    //                          heightMap,
    //                          cityModel,
    //                          mesh3D.DomainMarkers,
    //                          parameters.MeshResolution);

    // Generate height map function (used only for testing/visualization)
    //auto z = MeshSmoother::GenerateHeightMapFunction(_mesh2D, heightMap);

    // Generate mesh boundary (used only for testing/visualization)
    //dolfin::BoundaryMesh _boundary3D(_mesh3D, "exterior");

    // Write to files¨
    std::cout << "vc-generate-mesh: Writing to files..." << std::endl;
    dolfin::File(fileNamePrefix + "Mesh2D.pvd") << _mesh2D;
    // FIXME: Testing
    //dolfin::File(fileNamePrefix + "Mesh3D.pvd") << _mesh3D;
    //dolfin::File(fileNamePrefix + "Boundary.pvd") << _boundary3D;
    //dolfin::File(fileNamePrefix + "HeightMap.pvd") << *z;

    return 0;
}
