// vc-assimp-hello-world
// Orfeas Eleftheriou 2020

#include <fstream> //ifstream
#include <iostream>
#include <string>
#include <vector>

#include "Building.h"
#include "CityModel.h"
#include "JSON.h"
#include "Vector.h"

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

aiVector3D *generatePlaneVertices()
{
  // std::vector<aiVector3D> vertices;
  // vertices.push_back(aiVector3D(-1,-1,0));
  // vertices.push_back(aiVector3D(1,-1,0));
  // vertices.push_back(aiVector3D(1,1,0));
  // vertices.push_back(aiVector3D(-1,1,0));
  // return vertices;
  aiVector3D *vertices = new aiVector3D[4];
  vertices[0] = aiVector3D(-1, -1, 0);
  vertices[1] = aiVector3D(1, -1, 0);
  vertices[2] = aiVector3D(1, 1, 0);
  vertices[3] = aiVector3D(-1, 1, 0);
  return vertices;
}

aiVector3D *generatePlaneNormals()
{
  aiVector3D *normals = new aiVector3D[4];
  normals[0] = aiVector3D(0, 0, 1);
  normals[1] = aiVector3D(0, 0, 1);
  normals[2] = aiVector3D(0, 0, 1);
  normals[3] = aiVector3D(0, 0, 1);
  return normals;
}

aiVector3D *generatePlaneTangents()
{
  aiVector3D *tangents = new aiVector3D[4];
  tangents[0] = aiVector3D(1, 0, 0);
  tangents[1] = aiVector3D(1, 0, 0);
  tangents[2] = aiVector3D(1, 0, 0);
  tangents[3] = aiVector3D(1, 0, 0);
  return tangents;
}

aiVector3D *generatePlaneBitangents()
{
  aiVector3D *biTangents = new aiVector3D[4];
  biTangents[0] = aiVector3D(0, -1, 0);
  biTangents[1] = aiVector3D(0, -1, 0);
  biTangents[2] = aiVector3D(0, -1, 0);
  biTangents[3] = aiVector3D(0, -1, 0);
  return biTangents;
}

aiFace *generatePlaneFaces()
{
  aiFace *faces = new aiFace[2];

  // Create 1st face
  faces[0].mNumIndices = 3;
  faces[0].mIndices = new uint[3];

  faces[0].mIndices[0] = 0;
  faces[0].mIndices[1] = 1;
  faces[0].mIndices[2] = 2;

  // Create 2nd face
  faces[1].mNumIndices = 3;
  faces[1].mIndices = new uint[3];

  faces[1].mIndices[0] = 0;
  faces[1].mIndices[1] = 2;
  faces[1].mIndices[2] = 3;

  return faces;
}

aiMesh *generatePlaneMesh()
{
  aiMesh *mesh = new aiMesh();

  // Vertices
  mesh->mPrimitiveTypes = 4;
  mesh->mNumVertices = 4;
  mesh->mVertices = generatePlaneVertices();

  mesh->mName = "Plane";

  // Faces
  mesh->mNumFaces = 2;
  mesh->mFaces = generatePlaneFaces();

  // Normals
  mesh->mNormals = generatePlaneNormals();

  // Tangents
  mesh->mTangents = generatePlaneTangents();

  // BiTangents
  mesh->mBitangents = generatePlaneBitangents();

  return mesh;
}

aiMaterial *generatePlaneMaterial()
{
  aiMaterial *material = new aiMaterial();

  material->AddProperty(new aiVector3D(1, 1, 1), 1, AI_MATKEY_COLOR_DIFFUSE);
  return material;
}

aiScene *generatePlane()
{
  aiScene *plane = new aiScene();
  if (plane)
  {
    aiNode *rootNode = new aiNode();

    if (rootNode)
    {
      rootNode->mNumChildren = 1;
      rootNode->mChildren = new aiNode *[1];
      rootNode->mName = "RootNode";
      aiNode *subNode = new aiNode();

      if (subNode)
      {
        subNode->mParent = rootNode;
        subNode->mName = "Plane";
        subNode->mChildren = 0;
        subNode->mNumMeshes = 1;
        subNode->mMeshes = new uint[1];
        subNode->mMeshes[0] = 0;

        subNode->mTransformation = aiMatrix4x4(100, 0, 0, 0, 0, -1, 100, 0, 0,
                                               -100, -1, 0, 0, 0, 0, 1);
      }
      rootNode->mTransformation =
          aiMatrix4x4(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
      rootNode->mChildren[0] = subNode;
    }
    plane->mRootNode = rootNode;

    // Create a single solid material for the plane
    plane->mMaterials = new aiMaterial *[1];
    plane->mMaterials[0] = generatePlaneMaterial();
    plane->mNumMaterials = 1;

    // Store mesh
    plane->mFlags = 8;
    plane->mMeshes = new aiMesh *[1];
    plane->mNumMeshes = 1;
    plane->mMeshes[0] = generatePlaneMesh();
  }
  return plane;
}

int main(int argc, char *argv[])
{

  // Create an instance of the Importer class
  Assimp::Importer importer;
  // And have it read the given file with some example postprocessing
  // Usually - if speed is not the most important aspect for you - you'll
  // probably to request more postprocessing than we do in this example.

  // TODO: CHANGE THIS PATH
  const aiScene *cubeScene = importer.ReadFile(
      "cube.obj", aiProcess_CalcTangentSpace | aiProcess_Triangulate |
                      aiProcess_JoinIdenticalVertices | aiProcess_SortByPType);

  if (cubeScene)
  {
    std::cout << " valid cube scene with # " << cubeScene->mNumMeshes
              << " of meshes!" << std::endl;
  }
  else
  {
    std::cout << importer.GetErrorString() << std::endl;
  }

  std::cout << cubeScene->mRootNode->mNumChildren << std::endl;
  std::cout << cubeScene->mRootNode->mNumMeshes << std::endl;
  aiMesh *newMesh = cubeScene->mRootNode->mMeshes[0];

  /*
      Assimp::Exporter* TestExporter=new Assimp::Exporter();

      aiScene* planeScene = generatePlane();
      std::string fileName = "generatedmesh.obj";
      int result = TestExporter->Export(planeScene,"obj",fileName);
      if(result==0)
      {
          std::cout<<"File exported successfully: "<<fileName<<std::endl;
      }
      else
      {
          std::cout<<"There was an error in file export"<<std::endl;
          std::cout<<"Error:"<<TestExporter->GetErrorString()<<std::endl;
      }

  */

  return 0;
}
