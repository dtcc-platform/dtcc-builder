// vc-info
// Anders Logg 2019

#include <iostream>

#include "CityModel.h"
#include "CommandLine.h"
#include "HeightMap.h"
#include "JSON.h"
#include "LAS.h"
#include "Parameters.h"
#include "PointCloud.h"

using namespace DTCC;

void help() { std::cerr << "Usage: vc-info Data.[json,las]" << std::endl; }

int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (argc != 2)
  {
    help();
    return 1;
  }

  // Get filename
  const std::string fileName(argv[1]);

  // Check file type
  if (CommandLine::EndsWith(fileName, "json"))
  {
    const std::string type = JSON::ReadType(fileName);
    if (type == "Parameters")
    {
      Parameters parameters;
      JSON::Read(parameters, fileName);
      std::cout << parameters << std::endl;
    }
    else if (type == "HeightMap")
    {
      HeightMap heightMap;
      JSON::Read(heightMap, fileName);
      std::cout << heightMap << std::endl;
    }
    else if (type == "CityModel")
    {
      CityModel cityModel;
      JSON::Read(cityModel, fileName);
      std::cout << cityModel << std::endl;
    }
    else
    {
      std::cerr << "Unknown JSON type: \"" << type << "\"" << std::endl;
    }
  }
  else if (CommandLine::EndsWith(fileName, "las"))
  {
    PointCloud pointCloud;
    LAS::Read(pointCloud, fileName);
    std::cout << pointCloud << std::endl;
  }
  else
  {
    std::cerr << "Unhandled file type." << std::endl;
  }

  return 0;
}
