// Copyright (C) 2020 Anton J Olsson
// Licensed under the MIT LicenseS

#include "CommandLine.h"
#include "Logging.h"
#include "RoadNetworkGenerator.h"
#include <JSON.h>
#include <RoadNetwork.h>
#include <iostream>

using namespace DTCC;

void Help()
{
  std::cerr << "Usage: dtcc-generate-roadnetwork fileName.shp" << std::endl;
}

std::string GetDataDirectory(const std::string &filename);

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

  // Get RoadNetwork
  RoadNetwork network = RoadNetworkGenerator::GetRoadNetwork(shpFilename);

  // Serialize and write JSON file
  std::string dataDirectory = GetDataDirectory(shpFilename);
  json jsonNetwork;
  JSON::Serialize(network, jsonNetwork);
  JSON::Write(jsonNetwork, dataDirectory + "RoadNetwork.json");

  return 0;
}

// Get data directory given name of file located inside
std::string GetDataDirectory(const std::string &filename)
{
  int endPos = filename.find_last_of('/');
  return filename.substr(0, endPos + 1);
}
