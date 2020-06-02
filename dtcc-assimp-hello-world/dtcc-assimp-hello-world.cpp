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
#include <assimp/matrix4x4.h>


// Tasks:
// 1) Read file in str
// 2) Parse file's contents
// 3) Create data for assimp
// 4) Export with assimp

using namespace DTCC;

std::vector<aiVector3D> generateVertices()
{
    //These are relative locations to the placed Actor in the world
	// Vertices.Add(FVector(0, -100, 0)); //lower left - 0
	// Vertices.Add(FVector(0, -100, 100)); //upper left - 1
	// Vertices.Add(FVector(0, 100, 0)); //lower right - 2 
	// Vertices.Add(FVector(0, 100, 100)); //upper right - 3
	
	// Vertices.Add(FVector(100, -100, 0)); //lower front left - 4
	// Vertices.Add(FVector(100, -100, 100)); //upper front left - 5
 
	// Vertices.Add(FVector(100, 100, 100)); //upper front right - 6
	// Vertices.Add(FVector(100, 100, 0)); //lower front right - 7


    //Generate a triangle for now
    std::vector<aiVector3D> vertices;
    vertices.push_back(aiVector3D(0,-100,0));
    vertices.push_back(aiVector3D(0,-100,100));
    //vertices.push_back(aiVector3D(0,100,0));
    vertices.push_back(aiVector3D(0,100,100));

    //vertices.push_back(aiVector3D(100,-100,0));
    //vertices.push_back(aiVector3D(1000,-100,100));
//
    //vertices.push_back(aiVector3D(100,100,100));
    //vertices.push_back(aiVector3D(100,100,0));

    return vertices;
}

std::vector<aiFace> generateTriangles()
{
    std::vector<aiFace> triangles;
    aiFace face=aiFace();
    face.mNumIndices=3;
    face.mIndices=new uint[3];
    face.mIndices[0]=0;
    face.mIndices[1]=1;
    face.mIndices[2]=2;

    triangles.push_back(face);

    return triangles;
}

aiScene* generateCube()
{
    aiScene* scene=new aiScene();

    if(scene)
    {
        //Creating mesh
        aiMesh* mesh=new aiMesh();
        mesh->mPrimitiveTypes=4;
        const std::vector<aiVector3D> vertices=generateVertices();
        mesh->mNumVertices=vertices.size();

        mesh->mVertices=new aiVector3D[vertices.size()];
        mesh->mNormals=new aiVector3D[vertices.size()];
        mesh->mTangents=new aiVector3D[vertices.size()];
        

        //add vertices to actual mesh
        for(auto it=vertices.begin();it!=vertices.end();++it)
        {
            const aiVector3D& vertex = *it;

            mesh->mVertices[it-vertices.begin()]=aiVector3D(vertex.x,vertex.y,vertex.z);
            
            //std::string vertexstr = "x:"<<vertex.x<<" - y:"<<vertex.y<<" - z:"<<vertex.z;
            std::cout<<"added vertex to mesh:"<<"x:"<<vertex.x<<" - y:"<<vertex.y<<" - z:"<<vertex.z<<std::endl;
        }

        //Generate normals & tangents
        mesh->mNormals[0]=aiVector3D(0,0,1);
        mesh->mNormals[1]=aiVector3D(0,0,1);
        mesh->mNormals[2]=aiVector3D(0,0,1);

        mesh->mTangents[0]=aiVector3D(-1,0,0);
        mesh->mTangents[0]=aiVector3D(-1,0,0);
        mesh->mTangents[0]=aiVector3D(-1,0,0);

        //add triangles
        const std::vector<aiFace> triangles = generateTriangles();
        mesh->mFaces=new aiFace[triangles.size()];
        mesh->mNumFaces=triangles.size();

        for(auto it=triangles.begin();it!=triangles.end();++it)
        {
            aiFace& face=mesh->mFaces[it-triangles.begin()];

            face.mIndices=new uint[3];
            face.mNumIndices=3;
            face.mIndices[0] = 0;
            face.mIndices[1] = 1;
            face.mIndices[2] = 2;
        }

        //Create a new root node to store the mesh
        //scene->mRootNode=new aiNode();
        // scene->mRootNode->mMeshes=new uint[1];
        // scene->mRootNode->mMeshes[0]=0;
        // scene->mRootNode->mNumMeshes=1;

        aiNode* rootNode=new aiNode();

        if(rootNode)
        {
            rootNode->mNumChildren=1;
            rootNode->mChildren=new aiNode*[1]; //just 1 child for now
            aiNode* subNode=new aiNode();

            if(subNode)
            {
                subNode->mParent=rootNode;
                subNode->mChildren=0;
                subNode->mNumMeshes=1;
                subNode->mMeshes=new uint[1];
                subNode->mMeshes[0]=0;
            }
            rootNode->mChildren[0]=subNode;
        }
        scene->mRootNode=rootNode;

        //Store mesh
        scene->mFlags=8;
        scene->mMeshes=new aiMesh*[1];
        scene->mNumMeshes=1;
        scene->mMeshes[0]=mesh;
        
    }
    return scene;
}

aiVector3D* generatePlaneVertices()
{
    // std::vector<aiVector3D> vertices;
    // vertices.push_back(aiVector3D(-1,-1,0));
    // vertices.push_back(aiVector3D(1,-1,0));    
    // vertices.push_back(aiVector3D(1,1,0));    
    // vertices.push_back(aiVector3D(-1,1,0));
    // return vertices;
    aiVector3D* vertices=new aiVector3D[4];
    vertices[0]=aiVector3D(-1,-1,0);
    vertices[1]=aiVector3D(1,-1,0);    
    vertices[2]=aiVector3D(1,1,0);    
    vertices[3]=aiVector3D(-1,1,0);
    return vertices;
}

aiVector3D* generatePlaneNormals()
{
    aiVector3D* normals=new aiVector3D[4];
    normals[0]=aiVector3D(0,0,1);
    normals[1]=aiVector3D(0,0,1);
    normals[2]=aiVector3D(0,0,1);
    normals[3]=aiVector3D(0,0,1);
    return normals;
}

aiVector3D* generatePlaneTangents()
{
    aiVector3D* tangents=new aiVector3D[4];
    tangents[0]=aiVector3D(1,0,0);
    tangents[1]=aiVector3D(1,0,0);
    tangents[2]=aiVector3D(1,0,0);
    tangents[3]=aiVector3D(1,0,0);
    return tangents;
}

aiVector3D* generatePlaneBitangents()
{
    aiVector3D* biTangents=new aiVector3D[4];
    biTangents[0]=aiVector3D(0,-1,0);
    biTangents[1]=aiVector3D(0,-1,0);
    biTangents[2]=aiVector3D(0,-1,0);
    biTangents[3]=aiVector3D(0,-1,0);
    return biTangents;
}

aiFace* generatePlaneFaces()
{
    aiFace* faces=new aiFace[2];

    //Create 1st face
    faces[0].mNumIndices=3;
    faces[0].mIndices=new uint[3];

    faces[0].mIndices[0]=0;
    faces[0].mIndices[1]=1;
    faces[0].mIndices[2]=2;

    //Create 2nd face
    faces[1].mNumIndices=3;
    faces[1].mIndices=new uint[3];

    faces[1].mIndices[0]=0;
    faces[1].mIndices[1]=2;
    faces[1].mIndices[2]=3;

    return faces;
}

aiMesh* generatePlaneMesh()
{
    aiMesh* mesh=new aiMesh();

    //Vertices
    mesh->mPrimitiveTypes=4;
    mesh->mNumVertices=4;
    mesh->mVertices=generatePlaneVertices();

    //Faces
    mesh->mNumFaces=2;
    mesh->mFaces=generatePlaneFaces();

    //Normals
    mesh->mNormals=generatePlaneNormals();

    //Tangents
    mesh->mTangents=generatePlaneTangents();

    //BiTangents
    mesh->mBitangents=generatePlaneBitangents();

    return mesh;
}

aiScene* generatePlane()
{
    aiScene* plane=new aiScene();
    if(plane)
    {
        aiNode* rootNode=new aiNode();

        if(rootNode)
        {
            rootNode->mNumChildren=1;
            rootNode->mChildren=new aiNode*[1];
            aiNode* subNode=new aiNode();

            if(subNode)
            {
                subNode->mParent=rootNode;
                subNode->mName="Plane";
                subNode->mChildren=0;
                subNode->mNumMeshes=1;
                subNode->mMeshes=new uint[1];
                subNode->mMeshes[0]=0;

                subNode->mTransformation=aiMatrix4x4(100,0,0,0,0,-1,100,0,0,-100,-1,0,0,0,0,1);
            }
            rootNode->mTransformation=aiMatrix4x4(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);
            rootNode->mChildren[0]=subNode;
        }
        plane->mRootNode=rootNode;
        //Store mesh
        plane->mFlags=8;
        plane->mMeshes=new aiMesh*[1];
        plane->mNumMeshes=1;
        plane->mMeshes[0]=generatePlaneMesh();
    }
    return plane;
}

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


  //std::cout <<"----- json iteration -----"<<std::endl;

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

  //CityModel model;
  //JSON::Read(model,fileName);

  //bool firstBuildingRetrieved=false;
  //Building firstBuilding=Building();

  //for (auto it = model.Buildings.begin(); it != model.Buildings.end(); ++it)
  //{
  //    std::vector<Point2D> Points = it->Footprint.Points;
  //    for (auto p = Points.begin(); p != Points.end(); ++p)
  //    {
  //        std::cout<<"Point: "<<*p<<std::endl;
  //        if (!firstBuildingRetrieved)
  //        {
  //            firstBuilding.Footprint.Points.push_back(*p);
  //        }
  //    }
  //    //do-once flag
  //    if (!firstBuildingRetrieved)
  //    {
  //        firstBuilding.Height=it->Height;
  //        firstBuildingRetrieved=true;
  //    }
  //    std::cout<<"Moving on to next building!"<<std::endl;
  //}

  
  //Bottom-up design?
  //std::vector<aiVector3D> vertices;
  //std::cout<<"first building:"<<std::endl;
  //for (auto it = firstBuilding.Footprint.Points.begin(); it != firstBuilding.Footprint.Points.end(); ++it)
  //{
  //    std::cout<<"Point:"<<(*it)<<std::endl;
  //    aiVector3D newVertex;
  //    //newVertex.ai_real=0;
  //    newVertex.x=(*it).x;
  //    newVertex.y=(*it).y;
  //    newVertex.z=0; //no height yet?
  //    vertices.push_back(newVertex);
  //}

  aiMesh* mesh = new aiMesh();
  if(mesh)
  {
      /*mesh->mNumVertices=vertices.size();

      aiVector3D* meshVertices;
      meshVertices=(aiVector3D*)malloc(sizeof(aiVector3D) * mesh->mNumVertices);
      aiVector3D* firstElement = vertices.data();
      for (uint i = 0; i < mesh->mNumVertices; ++i)
      {
          meshVertices[i]=*firstElement;
          firstElement++;
      }*/

      mesh->mNumVertices=3;
      aiVector3D* meshVertices = new aiVector3D[3];
      meshVertices[0]=aiVector3D(0, 0, 0);
      meshVertices[1]=aiVector3D(0, 100, 0);
      meshVertices[2]=aiVector3D(0, 0, 100);

      mesh->mVertices=meshVertices;

    //   std::cout<< "Printing mesh vertices"<<std::endl;
    //   for (uint i = 0; i < mesh->mNumVertices; i++)
    //   {
    //       std::cout<<"X:"<< mesh->mVertices[i].x <<" - Y:" << mesh->mVertices[i].y << " - Z:" << mesh->mVertices[i].z << std::endl;
    //   }

      mesh->mNumFaces=1;

      mesh->mPrimitiveTypes=4;

	  //! Number of indices defining this face.
	  //! The maximum value for this member is #AI_MAX_FACE_INDICES.
	  //unsigned int mNumIndices;

	  //! Pointer to the indices array. Size of the array is given in numIndices.
	  //unsigned int* mIndices;
      /** The faces the mesh is constructed from.
       * Each face refers to a number of vertices by their indices.
       * This array is always present in a mesh, its size is given
       * in mNumFaces. If the #AI_SCENE_FLAGS_NON_VERBOSE_FORMAT
       * is NOT set each face references an unique set of vertices.
       */
      //C_STRUCT aiFace* mFaces;

      mesh->mFaces=new aiFace[1];
      mesh->mFaces[0].mNumIndices=3;
      mesh->mFaces[0].mIndices=new uint[3];
      mesh->mFaces[0].mIndices[0] = 0;
      mesh->mFaces[0].mIndices[1] = 1;
      mesh->mFaces[0].mIndices[2] = 2;
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


  // Create an instance of the Importer class
  Assimp::Importer importer;
  // And have it read the given file with some example postprocessing
  // Usually - if speed is not the most important aspect for you - you'll
  // probably to request more postprocessing than we do in this example.

  //TODO: CHANGE THIS PATH
//    const aiScene* cubeScene = importer.ReadFile("blender_plane.fbx",
//  	  aiProcess_CalcTangentSpace |
//  	  aiProcess_Triangulate |
//  	  aiProcess_JoinIdenticalVertices |
//  	  aiProcess_SortByPType);

//   if (cubeScene)
//   {
//       std::cout << " valid cube scene with # "<<cubeScene->mNumMeshes<<" of meshes!"<<std::endl;
//   }
//   else
//   {
//       std::cout <<importer.GetErrorString()<<std::endl;
//   }

//    Assimp::Exporter* Exporter=new Assimp::Exporter();
//    if (Exporter)
//    {
//        if (cubeScene)
//        {
//            std::cout<<"Valid scene"<<std::endl;
//            std::cout<<Exporter->Export(cubeScene, "fbx", "plane_export.fbx")<<std::endl;
//            std::cout << Exporter->Export(cubeScene, "assjson", "plane_export.assjson") << std::endl;
//            //std::cout << Exporter->Export(scene, "stl", "D:\\assimpexport.stl") << std::endl;
//            //std::cout << Exporter->Export(scene, "obj", "C:/GitHub_Projects/core/bin/assimpexport.obj") << std::endl;

// //           //std::cout << Exporter->ExportToBlob(scene,"stl") << std::endl;
// //           //std::cout << Exporter->Export(scene, "fbx", "D:/assimpexport.stl") << std::endl;
// //           std::cout << Exporter->Export(cubeScene, "fbx", "assimpexport.FBX") << std::endl;
// //           std::cout<<Exporter->Export(cubeScene,"assjson","asscube.assjson")<<std::endl;
// //           //std::cout<<Exporter->Export(scene,"assjson","assexport.assjson")<<std::endl;
// //           //std::cout<<Exporter->Export(scene,"fbx","testexp.fbx")<<std::endl;
// //           //std::cout<<"Error :"<<Exporter->GetErrorString()<<std::endl;

// //           //size_t ExportFormatCount = Exporter->GetExportFormatCount();
// //           //for (size_t i = 0; i < ExportFormatCount; i++)
// //           //{
// //           //    std::cout<<(Exporter->GetExportFormatDescription(i))->description<<std::endl;
// //           //}
//        }    
//    }
Assimp::Exporter* TestExporter=new Assimp::Exporter();
//aiScene* dummyScene = generateCube();
   
//    if(TestExporter)
//    {
//std::cout<<TestExporter->Export(dummyScene,"assjson","dummyscene.assjson")<<std::endl;
////std::cout<<TestExporter->Export(dummyScene,"fbx","dummyscene.fbx")<<std::endl;
//std::cout<<TestExporter->Export(dummyScene,"fbx","dummyscene.fbx")<<std::endl;
//std::cout<<"Error :"<<TestExporter->GetErrorString()<<std::endl;
   //}
 
    aiScene* planeScene = generatePlane();
    std::cout<<TestExporter->Export(planeScene,"assjson","generatedmesh.assjson")<<std::endl;
    std::cout<<"Error:"<<TestExporter->GetErrorString()<<std::endl;
    std::cout<<TestExporter->Export(planeScene,"fbx","generatedmesh.fbx")<<std::endl;
    //TestExporter->ExportToBlob(planeScene,"obj",1,nullptr);
    //std::cout<<TestExporter->Export(planeScene,"fbx","generatedmesh.fbx")<<std::endl;
    std::cout<<"Error:"<<TestExporter->GetErrorString()<<std::endl;

  return 0;
}
