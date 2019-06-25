// VirtualCity@Chalmers: vc-offset-geojson
// Anders Logg 2019
// Vasilis Naserentin 2019
// Offset in EPSG:3006 used in UE -148000 -6398600 (meters)
#include "GeoJSON.h"
#include "JSON.h"
#include <iostream>
#include <netcdf>
#include <string>
using namespace VirtualCity;

void Help()
{
  std::cerr << "Usage: vc-offset-geojson filein.geojson fileout.geojson x y "
            << std::endl;
}

/*
int detectmaxlevel(const nlohmann::json &ijson) {
  int maxlevel = 1;
  int curlevel = 0;
  while (maxlevel >= curlevel) {
    std::cout << "Trying for curlevel " << curlevel << std::endl;
    for (int k = 0; k < ijson["features"].size(); k++) {
      if (ijson["features"][k]["geometry"]["coordinates"][curlevel].size() >
          0) {
        std::cout
            << "Found "
            << ijson["features"][k]["geometry"]["coordinates"][curlevel].size()
            << std::endl;
        std::cout << "curlevel " << curlevel << " is valid" << std::endl;

        maxlevel++;
        curlevel++;
        std::cout << "New curlevel is " << curlevel << " and new maxlevel is "
                  << maxlevel << std::endl;
        break;
      }
      curlevel++;
    }
  }
  return maxlevel;
}
*/
int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (argc != 5)
  {
    Help();
    return 1;
  }
  int x, y;
  try
  {
    x = std::stod(argv[3]);
    y = std::stod(argv[4]);
  }
  catch (...)
  {
    std::cout << "Error; x and y must be integers." << std::endl;
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
  nlohmann::json json, jsoni;
  in >> json;
  jsoni = json;
  double offsetx = x;
  double offsety = y;
  std::cout << "JSON stuff" << std::endl;
  int maxlevel = GeoJSON::detectmaxlevel(json);
  std::cout << "Max level found is " << maxlevel << std::endl;
  std::cout << "Offsetting" << std::endl;

  for (int i = 0; i < maxlevel; i++)
  {
    for (int k = 0; k < json["features"].size(); k++)
    {
      for (int j = 0;
           j < json["features"][k]["geometry"]["coordinates"][i].size(); j++)
      {
        std::cout << i << " "
                  << json["features"][k]["geometry"]["coordinates"][i][j][0]
                  << std::endl;
        double oldx = json["features"][k]["geometry"]["coordinates"][i][j][0];
        double oldy = json["features"][k]["geometry"]["coordinates"][i][j][1];
        double newx = oldx + offsetx;
        double newy = oldy + offsety;
        jsoni["features"][k]["geometry"]["coordinates"][i][j][0] = newx;
        jsoni["features"][k]["geometry"]["coordinates"][i][j][1] = newy;
      }
    }
  }
  std::ofstream o(fileout);
  o << jsoni << std::endl;
  return 0;
}
