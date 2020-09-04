#include "CSV.h"
#include "JSON.h"
#include "CityJSON.h"
#include <iostream>
using namespace std;
using namespace DTCC;
void Help() { std::cerr << "Usage: vc-sandbox" << std::endl; }
int main(int argc, char *argv[])
{
  if (argc > 4)
  {
    Help();
    return 1;
  }
  // std::string varname = argv[3];
  // std::string fileout = argv[2];
  // std::string filein = argv[1];
  // nlohmann::json json;
  // in.read_header(io::ignore_extra_column, "vendor", "size", "speed");
  // std::string vendor; int size; double speed;
  CSV csv;
  //csv.Read("test.csv", true);
  CityJSON cityobj;
  JSON::Read(cityobj,"test.json");
  std::cout<<cityobj.CityObjects[0]<<std::endl;
  //DTCC::CityJSON::Read("HI");
}