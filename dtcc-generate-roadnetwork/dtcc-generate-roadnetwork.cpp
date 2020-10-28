// Copyright (C) 2020 Anton J Olsson
// Licensed under the MIT LicenseS

#include "CommandLine.h"
#include "Logging.h"
#include <Polygon.h>
#include <SHP.h>
#include <iostream>
#include <json.hpp>
#include <map>

using namespace DTCC;

void Help()
{
  std::cerr << "Usage: vc-generate-roadnetwork fileName.shp" << std::endl;
}

int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (argc != 2)
  {
    Help();
    return 1;
  }

  // Get .shp filename
  std::string shpFilename = std::string(argv[1]);

  Info(shpFilename);

  // Read road network data
  std::vector<Polygon> roadNetwork;
  nlohmann::json attributes;
  SHP::Read(roadNetwork, shpFilename, &attributes);

  return 0;
}
