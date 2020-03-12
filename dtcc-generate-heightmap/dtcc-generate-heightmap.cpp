// vc-generate-heightmap
// Anders Logg 2019

#include <iostream>
#include <string>
#include <vector>

#include "CommandLine.h"
#include "HeightMap.h"
#include "HeightMapGenerator.h"
#include "JSON.h"
#include "LAS.h"
#include "Parameters.h"

using namespace DTCC;

void Help()
{
  std::cerr << "Usage: vc-generate-heightmap Parameters.json" << std::endl;
}

int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (argc != 2)
  {
    Help();
    return 1;
  }

  // Read parameters
  Parameters parameters;
  JSON::Read(parameters, argv[1]);
  std::cout << parameters << std::endl;

  // Get data directory (add trailing slash just in case)
  const std::string dataDirectory = parameters.DataDirectory + "/";

  // Read point cloud data
  PointCloud pointCloud;
  for (auto const &f : CommandLine::ListDirectory(dataDirectory))
  {
    if (CommandLine::EndsWith(f, ".las") || CommandLine::EndsWith(f, ".laz"))
    {
      LAS::Read(pointCloud, dataDirectory + f);
      std::cout << pointCloud << std::endl;
    }
  }

  // Set domain size
  double xMin{}, yMin{}, xMax{}, yMax{};
  if (parameters.AutoDomain)
  {
    std::cout << "Automatically determining domain size:" << std::endl;
    xMin = pointCloud.XMin - parameters.X0;
    yMin = pointCloud.YMin - parameters.Y0;
    xMax = pointCloud.XMax - parameters.X0;
    yMax = pointCloud.YMax - parameters.Y0;
    std::cout << "  XMin: " << pointCloud.XMin << " --> " << xMin << std::endl;
    std::cout << "  YMin: " << pointCloud.YMin << " --> " << yMin << std::endl;
    std::cout << "  XMax: " << pointCloud.XMax << " --> " << xMax << std::endl;
    std::cout << "  YMax: " << pointCloud.YMax << " --> " << yMax << std::endl;
  }
  else
  {
    xMin = parameters.XMin;
    yMin = parameters.YMin;
    xMax = parameters.XMax;
    yMax = parameters.YMax;
  }

  // Generate height map
  HeightMap heightMap;
  HeightMapGenerator::GenerateHeightMap(
      heightMap, pointCloud, parameters.X0, parameters.Y0,
      xMin, yMin, xMax, yMax,
      parameters.HeightMapResolution);
  std::cout << heightMap << std::endl;

  // Write height map data
  JSON::Write(heightMap, dataDirectory + "HeightMap.json");

  return 0;
}
