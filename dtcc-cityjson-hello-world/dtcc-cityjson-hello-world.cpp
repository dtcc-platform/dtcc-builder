// Copyright (C) 2020 Vasilis Naserentin and Orfeas Eleftheriou
// Licensed under the MIT License

#include <iostream>
#include <string>
#include <vector>

#include "CommandLine.h"
#include "HeightMapGenerator.h"
#include "JSON.h"
#include "LAS.h"
#include "Logging.h"
#include "Parameters.h"
#include "cityjson/CityJSON.h"

using namespace DTCC;

void Help()
{
  Error("Usage: dtcc-generate-heightmap Parameters.json");
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
  /*Parameters parameters;
  JSON::Read(parameters, argv[1]);
  Info(parameters);*/

  nlohmann::json inputjson;

  CityJSON cityJson;
  JSON::Read(cityJson,argv[1]);

  // Debug Printing in case you want to access parsed data
  for(auto it:cityJson.CityObjects)
  {
    Progress(str(it));
  }

  for(auto it:cityJson.Vertices)
  {
    Progress(it.__str__());
  }

  //Testing Serialization to json
  JSON::Write(cityJson,"CityJsonTestExport.json");

  return 0;
}
