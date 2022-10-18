// Copyright (C) 2020 Orfeas Eleftheriou
// Licensed under the MIT License

#include <iostream>
#include <string>
#include <vector>

#include "CommandLine.h"
#include "JSON.h"
#include "LAS.h"
#include "Logging.h"
#include "Parameters.h"
#include "cityjson/CityJSON.h"

// Assimp includes
#include <assimp/Exporter.hpp>
#include <assimp/Importer.hpp>
#include <assimp/material.h>
#include <assimp/matrix4x4.h>
#include <assimp/mesh.h>
#include <assimp/postprocess.h>
#include <assimp/scene.h>
#include <assimp/vector3.h>

using namespace DTCC_BUILDER;

void Help()
{
  std::cerr << "Usage: dtcc-generate-cityjson-from-mesh Parameters.json"
            << std::endl;
}

std::vector<aiVector3D> getMeshVertices(const aiMesh &mesh)
{
  std::vector<aiVector3D> vertices;
  uint vertexNum = mesh.mNumVertices;

  for (uint i = 0; i < vertexNum; i++)
  {
    vertices.push_back(mesh.mVertices[i]);
  }
  return vertices;
}

std::vector<Point3D> getPointsFromVertices(const std::vector<aiVector3D> &vertices)
{
  std::vector<Point3D> points;
  for (auto &vertex : vertices)
  {
    points.push_back(Point3D(vertex.x, vertex.y, vertex.z));
  }
  return points;
}

std::vector<CityObject::Geometry::Boundary> getBoundariesFromMesh(const aiMesh &mesh)
{
  uint boundariesNum = mesh.mNumFaces;
  std::vector<CityObject::Geometry::Boundary> boundaries;

  for (uint i = 0; i < boundariesNum; i++)
  {
    uint faceIndices = mesh.mFaces[i].mNumIndices;
    CityObject::Geometry::Boundary boundary;
    for (uint j = 0; j < faceIndices; j++)
    {
      boundary.BoundariesIDs.push_back(mesh.mFaces[i].mIndices[j]);
    }
    boundaries.push_back(boundary);
  }

  return boundaries;
}

int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (argc != 2)
  {
    Help();
    return 1;
  }

  // Create an instance of the Importer class
  Assimp::Importer importer;
  // And have it read the given file with some example postprocessing
  // Usually - if speed is not the most important aspect for you - you'll
  // probably to request more postprocessing than we do in this example.
  std::string fileName = argv[1];
  // TODO: CHANGE THIS PATH
  const aiScene *importedScene = importer.ReadFile(
      fileName, aiProcess_CalcTangentSpace | aiProcess_Triangulate |
                    aiProcess_JoinIdenticalVertices | aiProcess_SortByPType);

  if (importedScene)
  {
    // std::cout << " valid cube scene with # " << cubeScene->mNumMeshes << " of
    // meshes!" << std::endl;

    // Get the first mesh
    std::cout << "Opening " << fileName << "..." << std::endl;

    aiMesh *mesh = importedScene->mMeshes[0];
    // Get vertices and print them on console
    std::vector<aiVector3D> vertices = getMeshVertices(*mesh);
    //printVertices(vertices);

    // Get the normals and print them on console
    //std::vector<aiVector3D> normals = getMeshNormals(*mesh);
    //printNormals(normals);

    std::vector<Point3D> points = getPointsFromVertices(vertices);
    std::vector<CityObject::Geometry::Boundary> boundaries =
        getBoundariesFromMesh(*mesh);

    CityObject cityObject;
    cityObject.ObjectType = CityObject::Building;
    cityObject.ObjectGeometry.Boundaries = boundaries;
    cityObject.ObjectGeometry.Type = CityObject::Geometry::Solid;

    CityJSON cityJson;
    cityJson.Vertices = points;
    cityJson.CityObjects.push_back(cityObject);
    fileName.append("_cityjson.json");
    JSON::Write(cityJson, fileName);
  }
  else
  {
    std::cout << importer.GetErrorString() << std::endl;
  }

  std::cout << importedScene->mRootNode->mNumChildren << std::endl;
  std::cout << importedScene->mRootNode->mNumMeshes << std::endl;

  return 0;
}
