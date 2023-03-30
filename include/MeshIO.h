// Copyright (C) 2020-2021 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_MESHIO_H
#define DTCC_MESHIO_H

// Assimp includes
#include <assimp/Exporter.hpp>
#include <assimp/Importer.hpp>
#include <assimp/material.h>
#include <assimp/mesh.h>
#include <assimp/scene.h>
#include <assimp/vector3.h>

#include "Logging.h"
#include "Mesh.h"
#include "Point.h"
#include "Simplex.h"
#include "Surface.h"

namespace DTCC
{

class MeshIO
{
public:
  static void Write(const Surface3D &surface,
                    std::string fileName,
                    std::string format,
                    bool YUp = false)
  {
    info("Assimp writing 3D surface to file " + fileName + " in format " +
         format);

    Assimp::Exporter *surfaceExporter = new Assimp::Exporter();
    aiScene *scene = new aiScene();
    if (scene)
    {
      aiNode *rootNode = new aiNode();
      rootNode->mNumMeshes = 1;
      rootNode->mMeshes = new uint[1];
      rootNode->mMeshes[0] = 0;
      scene->mRootNode = rootNode;

      scene->mMaterials = new aiMaterial *[1];
      scene->mMaterials[0] = defaultMaterial();
      scene->mNumMaterials = 1;

      scene->mFlags = 8;
      scene->mNumMeshes = 1;
      scene->mMeshes = new aiMesh *[1];
      scene->mMeshes[0] = buildAIMesh(surface.Vertices, surface.Faces, YUp);

      int result = surfaceExporter->Export(scene, format, fileName);
      if (result == 0)
      {
        info("File exported successfully: " + fileName);
      }
      else
      {
        error("There was an error in file export");
        error("Error:" + str(surfaceExporter->GetErrorString()));
      }
    }
    else
    {
      error("Assimp Error in creating scene");
    }
  }

private:
  static aiMaterial *defaultMaterial()
  {
    aiMaterial *material = new aiMaterial();

    material->AddProperty(new aiVector3D(1, 1, 1), 1, AI_MATKEY_COLOR_DIFFUSE);
    return material;
  }

  static aiMesh *buildAIMesh(std::vector<Point3D> vertices,
                             std::vector<Simplex2D> faces,
                             bool YUp)
  {
    aiMesh *mesh = new aiMesh();
    mesh->mPrimitiveTypes = aiPrimitiveType_TRIANGLE;
    mesh->mNumVertices = vertices.size();
    mesh->mNumFaces = faces.size();

    aiVector3D *aiVertices = new aiVector3D[mesh->mNumVertices];
    for (size_t i = 0; i < mesh->mNumVertices; i++)
    {
      auto p = vertices[i];
      if (YUp)
        aiVertices[i] = aiVector3D(p.x, p.z, p.y);
      else
        aiVertices[i] = aiVector3D(p.x, p.y, p.z);
    }
    mesh->mVertices = aiVertices;

    aiFace *aiFaces = new aiFace[mesh->mNumFaces];

    aiVector3D V1, V2, V3;
    double poly_sum;
    for (size_t i = 0; i < mesh->mNumFaces; i++)
    {
      aiFaces[i].mNumIndices = 3;
      aiFaces[i].mIndices = new uint[3];

      auto f = faces[i];
      // 30/3/2023 VN: Seems to produce inverted normals for walls and not for
      // ceilings, so we swapped V1 for V3 - Fix for JOSS review
      V3 = aiVertices[f.v0];
      V2 = aiVertices[f.v1];
      V1 = aiVertices[f.v2];

      // calculate triangle orientation using method from Geometry.h
      if (YUp)
      {
        poly_sum = (V2.x - V1.x) * (V2.z + V1.z) +
                   (V3.x - V2.x) * (V3.z + V2.z) +
                   (V1.x - V3.x) * (V1.z + V3.z);
      }
      else
      {
        poly_sum = (V2.x - V1.x) * (V2.y + V1.y) +
                   (V3.x - V2.x) * (V3.y + V2.y) +
                   (V1.x - V3.x) * (V1.y + V3.y);
      }

      // 30/3/2023 VN: Seems to produce inverted normals for walls and not for
      // ceilings, so we removed... - Fix for JOSS review if (poly_sum > 0)
      //{
      aiFaces[i].mIndices[0] = f.v0;
      aiFaces[i].mIndices[1] = f.v1;
      aiFaces[i].mIndices[2] = f.v2;
      // std::cout<<"1flipped"<<std::endl;
      /*   }
        else
        {
          aiFaces[i].mIndices[0] = f.v2;
          aiFaces[i].mIndices[1] = f.v1;
          aiFaces[i].mIndices[2] = f.v0;
          //std::cout<<"2not flipped"<<std::endl; */
      //}
    }
    mesh->mFaces = aiFaces;

    return mesh;
  }
};

} // namespace DTCC

#endif
