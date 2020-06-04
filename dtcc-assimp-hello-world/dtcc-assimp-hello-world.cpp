// vc-assimp-hello-world
// Orfeas Eleftheriou 2020

#include <fstream> //ifstream
#include <iostream>
#include <string>
#include <vector>

#include "Building.h"
#include "CityModel.h"
#include "JSON.h"
#include "Point.h"

#include <json.hpp>

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

aiVector3D* generatePlaneVertices()
{
  aiVector3D* vertices = new aiVector3D[4];
  vertices[0] = aiVector3D(-1, -1, 0);
  vertices[1] = aiVector3D(1, -1, 0);
  vertices[2] = aiVector3D(1, 1, 0);
  vertices[3] = aiVector3D(-1, 1, 0);
  return vertices;
}

aiVector3D* generatePlaneNormals()
{
  aiVector3D* normals = new aiVector3D[4];
  normals[0] = aiVector3D(0, 0, 1);
  normals[1] = aiVector3D(0, 0, 1);
  normals[2] = aiVector3D(0, 0, 1);
  normals[3] = aiVector3D(0, 0, 1);
  return normals;
}

aiVector3D* generatePlaneTangents()
{
  aiVector3D* tangents = new aiVector3D[4];
  tangents[0] = aiVector3D(1, 0, 0);
  tangents[1] = aiVector3D(1, 0, 0);
  tangents[2] = aiVector3D(1, 0, 0);
  tangents[3] = aiVector3D(1, 0, 0);
  return tangents;
}

aiVector3D* generatePlaneBitangents()
{
  aiVector3D* biTangents = new aiVector3D[4];
  biTangents[0] = aiVector3D(0, -1, 0);
  biTangents[1] = aiVector3D(0, -1, 0);
  biTangents[2] = aiVector3D(0, -1, 0);
  biTangents[3] = aiVector3D(0, -1, 0);
  return biTangents;
}

aiFace* generatePlaneFaces()
{
  aiFace* faces = new aiFace[2];

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

aiMesh* generatePlaneMesh()
{
  aiMesh* mesh = new aiMesh();

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

aiMaterial* generatePlaneMaterial()
{
  aiMaterial* material = new aiMaterial();

  material->AddProperty(new aiVector3D(1, 1, 1), 1, AI_MATKEY_COLOR_DIFFUSE);
  return material;
}

aiScene* generatePlane()
{
  aiScene *plane = new aiScene();
  if (plane)
  {
    aiNode *rootNode = new aiNode();

    if (rootNode)
    {
      rootNode->mNumChildren = 1;
      rootNode->mChildren = new aiNode*[1];
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
    plane->mMaterials = new aiMaterial*[1];
    plane->mMaterials[0] = generatePlaneMaterial();
    plane->mNumMaterials = 1;

    // Store mesh
    plane->mFlags = 8;
    plane->mMeshes = new aiMesh*[1];
    plane->mNumMeshes = 1;
    plane->mMeshes[0] = generatePlaneMesh();
  }
  return plane;
}

void ExportTestPlane() 
{
  Assimp::Exporter* TestExporter = new Assimp::Exporter();

  aiScene* planeScene = generatePlane();
  std::string fileName = "generatedmesh.obj";
  int result = TestExporter->Export(planeScene, "obj", fileName);
  if (result == 0)
  {
    std::cout << "File exported successfully: " << fileName << std::endl;
  }
  else
  {
    std::cout << "There was an error in file export" << std::endl;
    std::cout << "Error:" << TestExporter->GetErrorString() << std::endl;
  }
}

nlohmann::json Read(std::string fileName)
{
  //Info("JSON: Reading from file " + fileName + "...");
  std::ifstream f(fileName);
  if (!f)
  {
    std::cerr << "Unable to read from file"<<fileName<<std::endl;
  }
    
  nlohmann::json json{};
  f >> json;
  return json;
}

std::vector<std::vector<uint>> GetBoundariesFromGeometry(nlohmann::json geometryJson)
{
  std::vector<std::vector<uint>> boundaries;
  uint boundariesNum = geometryJson["boundaries"][0].size();

  for(uint index=0;index<boundariesNum;index++)
  {
    auto face = geometryJson["boundaries"][0][index];
    for(auto faceIt=face.begin();faceIt!=face.end();++faceIt)
    {
      //std::cout << "adding face #"<<index<<std::endl;
      std::vector<uint> faceIndices=faceIt.value();
      boundaries.push_back(faceIndices);
    }
  }
  
  return boundaries;
}

std::vector<aiVector3D> GetVerticesFromJson(nlohmann::json jsonFile) 
{
  std::vector<aiVector3D> vertices;
  uint num = jsonFile["vertices"].size();
  for(uint i=0;i<num;i++)
  {
    std::cout << jsonFile["vertices"][i] << std::endl;
    auto vertex = jsonFile["vertices"][i];
    std::vector<double> tempVertexArray;
    for(auto it=vertex.begin();it!=vertex.end();++it)
    {
        tempVertexArray.push_back((double)*it);
    }
    vertices.push_back(aiVector3D(tempVertexArray[0],tempVertexArray[1],tempVertexArray[2]));
  }
  return vertices;
}

aiScene* CreateMesh(const std::vector<aiVector3D>& vertices, const std::vector<std::vector<uint>>& faces)
{
    aiScene* scene= new aiScene();
    if (scene)
    {
        aiNode* rootNode=new aiNode();

        if (rootNode)
        {
            rootNode->mNumChildren=1;
            rootNode->mChildren = new aiNode *[1];
            rootNode->mName = "RootNode";
            aiNode* subNode=new aiNode();

            if (subNode)
            {
              subNode->mParent=rootNode;
              subNode->mName = "BuildingTest";
              subNode->mChildren=0;
              subNode->mNumMeshes=1;
              subNode->mMeshes = new uint[1];
              subNode->mMeshes[0]=0;
            }

            rootNode->mChildren[0]=subNode;
        }

        scene->mRootNode=rootNode;
        scene->mMaterials = new aiMaterial *[1];
        scene->mMaterials[0]=generatePlaneMaterial();
        scene->mNumMaterials=1;

        scene->mFlags=8;
        scene->mMeshes = new aiMesh*[1];
        scene->mNumMeshes=1;

        aiVector3D *meshVertices = new aiVector3D[vertices.size()];
        for (auto it = vertices.begin(); it != vertices.end(); ++it)
        {
          meshVertices[it - vertices.begin()] = *it;
        }

        aiFace *aiFaces = new aiFace[faces.size()];

        for (auto faceIt = faces.begin(); faceIt != faces.end(); ++faceIt)
        {
          (aiFaces[faceIt - faces.begin()]).mNumIndices = faceIt->size();
          (aiFaces[faceIt-faces.begin()]).mIndices=new uint(faceIt->size());
          for(auto it=faceIt->begin();it!=faceIt->end();++it)
          {
            (aiFaces[faceIt - faces.begin()]).mIndices[it-faceIt->begin()]=*it;
          }
        }
        aiMesh* mesh= new aiMesh();
        mesh->mName = "TestBuilding";
        //vertices
        mesh->mPrimitiveTypes=4;
        mesh->mNumVertices=vertices.size();
        mesh->mVertices=meshVertices;

        //faces
        mesh->mNumFaces=faces.size();
        mesh->mFaces=aiFaces;

        scene->mMeshes[0]=mesh;
    }
    return scene;
}

void ExportMesh(const std::vector<aiVector3D>& vertices, const std::vector<std::vector<uint>>& faces, const std::string& extension, const std::string& fileName)
{
  aiScene* sceneToExport = CreateMesh(vertices,faces);
  Assimp::Exporter *TestExporter = new Assimp::Exporter();
  int result = TestExporter->Export(sceneToExport, extension, fileName,
                                    aiProcess_ForceGenNormals);
  if (result == 0)
  {
    std::cout << "File exported successfully: " << fileName << std::endl;
  }
  else
  {
    std::cout << "There was an error in file export" << std::endl;
    std::cout << "Error:" << TestExporter->GetErrorString() << std::endl;
  }
}

int main(int argc, char *argv[])
{
  // Check command-line arguments
  // std::cout << "Hello World!"<<std::endl;

  nlohmann::json js=Read("building.json");
  //std::cout<<js["type"]<<std::endl;
  auto geometry = js["CityObjects"]["UUID_d281adfc-4901-0f52-540b-4cc1a9325f82"]["geometry"];
  std::vector<std::vector<uint>> boundaries;
  boundaries=GetBoundariesFromGeometry(geometry[0]);

  //debug print boundaries
  /*for(auto boundary=boundaries.begin();boundary!=boundaries.end();++boundary)
  {
    for(auto faceIt=boundary->begin();faceIt!=boundary->end();++faceIt)
    {
      std::cout<<*faceIt<<std::endl;
    }
  }*/
  
  std::cout << "num of faces:"<<boundaries.size()<<std::endl;
  std::cout<<"vertices:"<<std::endl<<js["vertices"]<<std::endl;

  std::vector<aiVector3D> vertices = GetVerticesFromJson(js);
  
  ExportMesh(vertices,boundaries,"fbx","TestBuilding_ForcedNormals.fbx");


  return 0;
}
