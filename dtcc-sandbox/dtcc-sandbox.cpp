// Copyright (C) 2020-2021 Anders Logg, Anton J Olsson
// Licensed under the MIT License

#include <string>
#include <vector>

// Needs to come before JSON (nlohmann) include because of sloppy
// namespacing in VTK (typedef detail)...

// DTCC includes
#include "JSON.h"
#include "Mesh.h"

using namespace DTCC_BUILDER;

void Help() { error("Usage: dtcc-generate-citymodel Parameters.json"); }


int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (argc != 2)
  {
    Help();
    return 1;
  }

  // Read parameters
  Parameters p;
  Mesh3D m;
  JSON::Read(m, argv[1]);
  info(m);
}
