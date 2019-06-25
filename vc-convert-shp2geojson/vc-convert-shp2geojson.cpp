// VirtualCity@Chalmers: vc-convert-shp2geojson
// Anders Logg 2019
// Vasilis Naserentin 2019
// Requires ogr2ogr located in /usr/bin/ in VC Docker image
// example ogr2ogr usage:
//./ogr2ogr -f GeoJSON  -s_srs EPSG:3006 -t_srs EPSG:3006
/// tmp/PropertyMap.geojson /tmp/PropertyMap.shp

#include "JSON.h"
#include <iostream>
#include <string>

void Help()
{
  std::cerr << "Usage: vc-convert-shp2geojson in.shp out.geojson" << std::endl;
}

int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (argc != 3)
  {
    Help();
    return 1;
  }
  std::string fileout = argv[2];

  std::string filein = argv[1];
  std::ifstream in(filein);
  if (in.fail())
  {
    std::cout << "Input file issue; check input file exists." << std::endl;
    return 1;
  }
  in.close();
  std::ifstream ogr("/usr/bin/ogr2ogr");
  if (ogr.fail())
  {
    std::cout << "ogr2ogr binary not found. Please install gdal-bin."
              << std::endl;
    return 1;
  }
  ogr.close();
  system(("ogr2ogr -f GeoJSON  -s_srs EPSG:3006 -t_srs EPSG:3006 " + fileout +
          " " + filein)
             .c_str());
  return 0;
}
