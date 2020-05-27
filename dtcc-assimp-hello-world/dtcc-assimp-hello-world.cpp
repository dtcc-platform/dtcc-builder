// vc-assimp-hello-world
// Orfeas Eleftheriou 2020 

#include <iostream>
#include <iostream>
#include <fstream> //ifstream
#include <string>
#include <vector>

#include "JSON.h"
#include "Building.h"
#include "CityModel.h"
#include "Point.h"

//#include <json.hpp>
#include <assimp/Importer.hpp>      
#include <assimp/Exporter.hpp>
#include <assimp/scene.h>           
#include <assimp/mesh.h>
#include <assimp/vector3.h>
#include <assimp/postprocess.h>     



// Tasks:
// 1) Read file in str
// 2) Parse file's contents
// 3) Create data for assimp
// 4) Export with assimp

using namespace DTCC;

int main(int argc, char *argv[])
{
  // Check command-line arguments
  std::cout << "Hello World!"<<std::endl;
  //Assimp::Importer Importer;
  std::string fileName = "CityModel.json";
  //std::ifstream infile("file.json");
  //std::ifstream inFile(fileName);
  std::ifstream inFile(fileName);
  
  //nlohmann::json json;
  //inFile >> json;

  //std::string serialized_str=json.dump();

  //Testing outputs
  //std::cout << "infile:" << inFile <<std::endl;
  //std::cout <<"nlohman json:"<<serialized_str<<std::endl;


  std::cout <<"----- json iteration -----"<<std::endl;

  //int maxPrints=10;
  //int currentPrints=0;
  //for (auto it = json.begin(); it != json.end() /*&& currentPrints<maxPrints*/; ++it)
  //{
  //    std::cout << "key:" << it.key() <<" : "<< it.value() <<std::endl;
  //    //currentPrints++;
  //}

  //std::vector<Building> buildings = json["Buildings"].get<std::vector<Building>>();

  //for (auto it = buildings.begin(); it != buildings.end(); ++it)
  //{
  //    //std::cout<<*it<<std::endl;
  //}

  CityModel model;
  JSON::Read(model,fileName);

  bool firstBuildingRetrieved=false;
  Building firstBuilding=Building();

  for (auto it = model.Buildings.begin(); it != model.Buildings.end(); ++it)
  {
      std::vector<Point2D> Points = it->Footprint.Points;
      for (auto p = Points.begin(); p != Points.end(); ++p)
      {
          std::cout<<"Point: "<<*p<<std::endl;
          if (!firstBuildingRetrieved)
          {
              firstBuilding.Footprint.Points.push_back(*p);
          }
      }
      //do-once flag
      if (!firstBuildingRetrieved)
      {
          firstBuilding.Height=it->Height;
          firstBuildingRetrieved=true;
      }
      std::cout<<"Moving on to next building!"<<std::endl;
  }

  
  //Bottom-up design?
  std::vector<aiVector3D> vertices;
  std::cout<<"first building:"<<std::endl;
  for (auto it = firstBuilding.Footprint.Points.begin(); it != firstBuilding.Footprint.Points.end(); ++it)
  {
      std::cout<<"Point:"<<(*it)<<std::endl;
      aiVector3D newVertex;
      //newVertex.ai_real=0;
      newVertex.x=(*it).x;
      newVertex.y=(*it).y;
      newVertex.z=0; //no height yet?
      vertices.push_back(newVertex);
  }

  aiMesh* mesh = new aiMesh();
  if(mesh)
  {
      mesh->mNumVertices=vertices.size();

      aiVector3D* meshVertices;
      meshVertices=(aiVector3D*)malloc(sizeof(aiVector3D) * mesh->mNumVertices);
      aiVector3D* firstElement = vertices.data();
      for (uint i = 0; i < mesh->mNumVertices; ++i)
      {
          meshVertices[i]=*firstElement;
          firstElement++;
      }

      mesh->mVertices=meshVertices;

      std::cout<< "Printing mesh vertices"<<std::endl;
      for (uint i = 0; i < mesh->mNumVertices; i++)
      {
          std::cout<<"X:"<< mesh->mVertices[i].x <<" - Y:" << mesh->mVertices[i].y << std::endl;
      }
  }

  aiNode* rootNode = new aiNode();
  if (rootNode)
  {
      rootNode->mNumMeshes=1; //for now
      //rootNode->mMeshes=1; //more info: https://github.com/assimp/assimp/blob/master/include/assimp/scene.h#L128
      rootNode->mMeshes=(uint*)malloc(sizeof(uint)*1);
      rootNode->mMeshes[0]=0;
  }

  aiScene* scene = new aiScene();
  if (scene)
  {
      scene->mNumMeshes=1; //for now
      scene->mRootNode=rootNode;
      scene->mMeshes = (aiMesh**)malloc(sizeof(aiMesh)*1); //1 mesh for now
      scene->mMeshes[0]=mesh;
      
  }

  /*size_t formatsNum= Assimp::Exporter::GetExportFormatCount();
  for (size_t i = 0; i < formatsNum; i++)
  {
      std::cout<<"file extension:"<<Assimp::Exporter::GetExportFormatDescription(i)->fileExtension()<<std::endl;
  }*/

  Assimp::Exporter* Exporter=new Assimp::Exporter();
  if (Exporter)
  {
      if (scene)
      {
          std::cout<<"Valid scene"<<std::endl;
          //std::cout<<Exporter->Export(scene, "fbx", "C:/GitHub_Projects/core/bin/assimpexport")<<std::endl;
          std::cout << Exporter->Export(scene, "stl", "C:/GitHub_Projects/core/bin/assimpexport") << std::endl;
          std::cout<<"Error :"<<Exporter->GetErrorString()<<std::endl;

          size_t ExportFormatCount = Exporter->GetExportFormatCount();
          for (size_t i = 0; i < ExportFormatCount; i++)
          {
              std::cout<<(Exporter->GetExportFormatDescription(i))->description<<std::endl;
          }
      }    
  }

 /* std::vector<XYZ>& points = m_pointMaps[0];

  aiScene scene;
  scene.mRootNode = new aiNode();

  scene.mMeshes = new aiMesh * [1];
  scene.mMeshes[0] = nullptr;
  scene.mNumMeshes = 1;

  scene.mMaterials = new aiMaterial * [1];
  scene.mMaterials[0] = nullptr;
  scene.mNumMaterials = 1;

  scene.mMaterials[0] = new aiMaterial();

  scene.mMeshes[0] = new aiMesh();
  scene.mMeshes[0]->mMaterialIndex = 0;

  scene.mRootNode->mMeshes = new unsigned int[1];
  scene.mRootNode->mMeshes[0] = 0;
  scene.mRootNode->mNumMeshes = 1;

  auto pMesh = scene.mMeshes[0];

  long numValidPoints = points.size() - sumCutPixels;

  pMesh->mVertices = new aiVector3D[numValidPoints];
  pMesh->mNumVertices = numValidPoints;

  int i = 0;
  for (XYZ& p : points)
  {
	  if (isnan(p.x) or isnan(p.y) or isnan(p.z))
		  continue;

	  pMesh->mVertices[i] = aiVector3D(p.x, p.y, p.z);
	  ++i;
  }*/

  return 0;
}
