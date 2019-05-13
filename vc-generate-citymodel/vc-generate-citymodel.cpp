// VirtualCity@Chalmers: vc-generate-citymodel
// Anders Logg 2019

#include <iostream>

#include "CommandLine.h"
#include "OSM.h"
#include "SHP.h"
#include "JSON.h"
#include "CityModel.h"
#include "CityModelGenerator.h"

using namespace VirtualCity;

void Help()
{
  std::cerr << "Usage: vc-generate-citymodel PropertyMap.[shp/osm]"
            << " HeightMap.json Parameters.json" << std::endl;
}

int main(int argc, char* argv[])
{
  // Check command-line arguments
  if (argc != 4)
  {
    Help();
    return 1;
  }

  // Get filenames
  const std::string fileNamePropertyMap(argv[1]);
  const std::string fileNameHeightMap(argv[2]);
  const std::string fileNameParameters(argv[3]);

  // FIXME: OSM untested and needs some further work
  if (CommandLine::EndsWith(fileNamePropertyMap, ".osm"))
  {
    std::cout << "Warning: OSM data untested and likely not working atm."
              << std::endl;
  }

  // Read polygons
  std::vector<Polygon> polygons;
  if (CommandLine::EndsWith(fileNamePropertyMap, ".shp"))
    SHP::Read(polygons, fileNamePropertyMap);
  else
    OSM::Read(polygons, fileNamePropertyMap);

  // Read height map
  HeightMap heightMap;
  JSON::Read(heightMap, fileNameHeightMap);
  std::cout << heightMap << std::endl;

  // Read parameters from file
  Parameters parameters;
  JSON::Read(parameters, fileNameParameters);
  std::cout << parameters << std::endl;

  // Generate city model
  CityModel cityModel;
  CityModelGenerator::GenerateCityModel(cityModel,
                                        polygons,
                                        heightMap,
                                        parameters.X0, parameters.Y0,
                                        parameters.XMin, parameters.YMin,
                                        parameters.XMax, parameters.YMax,
                                        parameters.MinimalBuildingDistance);
  std::cout << cityModel << std::endl;

  // Write to file
  JSON::Write(cityModel, "CityModel.json");

  return 0;
}
