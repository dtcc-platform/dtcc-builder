// vc-sandbox
// Just a sandbox code for users to experiment with different API calls for
// VCCore Vasilis Naserentin 2019
// Licensed under the MIT License

#include "CSV.h"
#include "JSON.h"
#include <iostream>

using namespace std;

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
  DTCC::CSV csv;
  csv.Read("test.csv", true);
}
