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
#include "CityJSON.h"

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


  return 0;
}
