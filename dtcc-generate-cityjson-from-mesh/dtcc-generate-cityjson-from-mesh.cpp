// Copyright (C) 2020 Orfeas Eleftheriou
// Licensed under the MIT License

#include <iostream>
#include <string>
#include <vector>

#include "CommandLine.h"
#include "HeightMapGenerator.h"
#include "JSON.h"
#include "LAS.h"
#include "Parameters.h"
#include "CityJSON.h"
#include "Logging.h"

// Assimp includes
#include <assimp/Exporter.hpp>
#include <assimp/Importer.hpp>
#include <assimp/material.h>
#include <assimp/matrix4x4.h>
#include <assimp/mesh.h>
#include <assimp/postprocess.h>
#include <assimp/scene.h>
#include <assimp/vector3.h>

using namespace DTCC;

void Help()
{
  std::cerr << "Usage: vc-generate-heightmap Parameters.json" << std::endl;
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
  Assimp::Importer importer;
  // TODO: CHANGE THIS PATH
  const aiScene *cubeScene = importer.ReadFile(
      argv[1], aiProcess_CalcTangentSpace | aiProcess_Triangulate |
                      aiProcess_JoinIdenticalVertices | aiProcess_SortByPType);

  if(!cubeScene)
  {
    std::cout << "Failed to import file:"<<importer.GetErrorString()<<std::endl;
  }

  //nlohmann::json inputjson;

  //CityJSON cityJson;
  //JSON::Read(cityJson,argv[1]);

  //// Debug Printing in case you want to access parsed data
  //for(auto it:cityJson.CityObjects)
  //{
  //  std::cout<<it<<std::endl;
  //}

  //for(auto it:cityJson.Vertices)
  //{
  //  std::cout<<it.__str__()<<std::endl;
  //}

  ////Testing Serialization to json
  //JSON::Write(cityJson,"CityJsonTestExport.json");

  return 0;
}

