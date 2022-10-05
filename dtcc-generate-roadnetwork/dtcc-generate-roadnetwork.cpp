// Copyright (C) 2020-2021 Anton J Olsson
// Licensed under the MIT LicenseS

#include "CommandLine.h"
#include "Logging.h"
#include "RoadNetworkGenerator.h"
#include <JSON.h>
#include <Polygon.h>
#include <SHP.h>
#include <datamodel/RoadNetwork.h>
#include <Polyfix.h>
#include <datamodel/RoadNetwork.h>
#include <iostream>

using namespace DTCC;

void Help()
{
  error("Usage: dtcc-generate-roadnetwork fileName.shp");
}

std::string GetDataDirectory(const std::string &filename);
void OffsetCoordinates(RoadNetwork &network,
                       const Point2D &origin,
                       Polygon &networkPolygon);
Polygon GetPolygon(const RoadNetwork &network);

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
  info(shpFilename);

  // Get RoadNetwork
  RoadNetwork network = RoadNetworkGenerator::GetRoadNetwork(shpFilename);

  // Offset coordinates
  Polygon networkPolygon = GetPolygon(network);
  Point2D origin = BoundingBox2D(networkPolygon).P;
  OffsetCoordinates(network, origin, networkPolygon);

  // Serialize and write JSON file
  std::string dataDirectory = GetDataDirectory(shpFilename);
  json jsonNetwork;
  JSON::Serialize(network, jsonNetwork, origin);
  JSON::Write(jsonNetwork, dataDirectory + "RoadNetwork.json");

  return 0;
}

Polygon GetPolygon(const RoadNetwork &network)
{
  Polygon networkPolygon;
  networkPolygon.Vertices = network.Vertices;
  return networkPolygon;
}

void OffsetCoordinates(RoadNetwork &network,
                       const Point2D &origin,
                       Polygon &networkPolygon)
{
  networkPolygon.SetOrigin(origin);
  network.Vertices = networkPolygon.Vertices;
}

// Get data directory given name of file located inside
std::string GetDataDirectory(const std::string &filename)
{
  int endPos = filename.find_last_of('/');
  return filename.substr(0, endPos + 1);
}
