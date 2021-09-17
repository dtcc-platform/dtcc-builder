// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#include <dolfin.h>

#include "CityModelGenerator.h"
#include "CommandLine.h"
#include "GridField.h"
#include "JSON.h"
#include "LaplacianSmoother.h"
#include "Logging.h"
#include "Mesh.h"
#include "MeshGenerator.h"
#include "MeshProcessor.h"
#include "OBJ.h"
#include "Parameters.h"
#include "VTK.h"
#include "VertexSmoother.h"

using namespace DTCC;

void Help() { Error("Usage: dtcc-generate-simulation-mesh Parameters.json"); }

// Generate surface meshes (non-matching, used for visualization)
void GenerateSurfaceMeshes(const CityModel &cityModel,
                           const GridField2D &dtm,
                           const Parameters &p)
{
  // Set data directory
  const std::string dataDirectory{p.DataDirectory + "/"};

  // Get origin (for serialization purposes)
  Point2D origin({p.X0, p.Y0});

  // Generate mesh for ground and buildings
  Surface3D groundSurface;
  std::vector<Surface3D> buildingSurfaces;
  MeshGenerator::GenerateSurfaces3D(groundSurface, buildingSurfaces, cityModel,
                                    dtm, p.MeshResolution, p.FlatGround);

  // Merge building surfaces
  Surface3D buildingSurface;
  MeshProcessor::MergeSurfaces3D(buildingSurface, buildingSurfaces);

  // Write to file
  if (p.WriteJSON)
  {
    JSON::Write(groundSurface, dataDirectory + "GroundSurface.json", origin);
    JSON::Write(buildingSurface, dataDirectory + "BuildingSurface.json",
                origin);
  }

  // Write data for debugging and visualization
  if (p.WriteVTK)
  {
    VTK::Write(groundSurface, dataDirectory + "GroundSurface.vtu");
    VTK::Write(buildingSurface, dataDirectory + "BuildingSurface.vtu");
  }
}

// Generate volume meshes (matching, used for simulation)
void GenerateVolumeMeshes(CityModel &cityModel,
                          const GridField2D &dtm,
                          const Parameters &p)
{
  // Set data directory
  const std::string dataDirectory{p.DataDirectory + "/"};

  // Get origin (for serialization purposes)
  Point2D origin({p.X0, p.Y0});

  // Step 1: Generate city model (and elevation model).
  // This step is handled by dtcc-generate-citymodel and
  // we assume that the data has already been generated.

  // Step 2.1: Merge building footprints
  {
    Timer timer("Step 2.1: Merge building footprints");
    CityModelGenerator::SimplifyCityModel(cityModel, p.MinBuildingDistance,
                                          p.MinVertexDistance);
    Info(cityModel);
  }

  // Step 2.2: Clean building footprints
  {
    Timer timer("Step 2.2: Clean building footprints");
    CityModelGenerator::CleanCityModel(cityModel, p.MinVertexDistance);
    Info(cityModel);
  }

  // Step 2.3: Recompute building heights from building points and DTM
  {
    Timer timer("Step 2.3: Clean building footprints");
    CityModelGenerator::ComputeBuildingHeights(
        cityModel, dtm, p.GroundPercentile, p.RoofPercentile);
    Info(cityModel);
  }

  // Write to file
  if (p.WriteJSON)
  {
    JSON::Write(cityModel, dataDirectory + "CityModelSimple.json", origin);
  }

  // Step 3.1: Generate 2D mesh
  Mesh2D mesh2D;
  {
    Timer timer("Step 3.1: Generate 2D mesh");
    MeshGenerator::GenerateMesh2D(mesh2D, cityModel, dtm.Grid.BoundingBox,
                                  p.MeshResolution);
    Info(mesh2D);
  }

  // Write data for debugging and visualization
  if (p.WriteVTK)
  {
    VTK::Write(mesh2D, dataDirectory + "Step31Mesh.vtu");
  }

  // Step 3.2: Generate 3D mesh (layer 2D mesh)
  size_t numLayers{};
  Mesh3D mesh;
  {
    Timer timer("Step 3.2: Generate 3D mesh (layer 2D mesh)");
    numLayers = MeshGenerator::GenerateMesh3D(mesh, mesh2D, p.DomainHeight,
                                              p.MeshResolution);
    Info(mesh);
  }

  // Write data for debugging and visualization
  if (p.WriteVTK)
  {
    Surface3D boundary;
    MeshProcessor::ExtractBoundary3D(boundary, mesh);
    VTK::Write(mesh, dataDirectory + "Step32Mesh.vtu");
    VTK::Write(boundary, dataDirectory + "Step32Boundary.vtu");
  };

  // Step 3.3: Smooth 3D mesh (set ground height)
  double topHeight{};
  {
    Timer timer("Step 3.3: Smooth 3D mesh (set ground height)");
    topHeight = dtm.Mean() + p.DomainHeight;
    LaplacianSmoother::SmoothMesh3D(mesh, cityModel, dtm, topHeight, false);
    Info(mesh);
  }

  // Write data for debugging and visualization
  if (p.WriteVTK)
  {
    Surface3D boundary;
    MeshProcessor::ExtractBoundary3D(boundary, mesh);
    VTK::Write(mesh, dataDirectory + "Step33Mesh.vtu");
    VTK::Write(boundary, dataDirectory + "Step33Boundary.vtu");
  }

  // Step 3.4: Trim 3D mesh (remove building interiors)
  {
    Timer timer("Step 3.4: Trim 3D mesh (remove building interiors)");
    MeshGenerator::TrimMesh3D(mesh, mesh2D, cityModel, numLayers);
    Info(mesh);
  }

  // Write data for debugging and visualization
  if (p.WriteVTK)
  {
    Surface3D boundary;
    MeshProcessor::ExtractBoundary3D(boundary, mesh);
    VTK::Write(mesh, dataDirectory + "Step34Mesh.vtu");
    VTK::Write(boundary, dataDirectory + "Step34Boundary.vtu");
  }

  // Step 3.5: Smooth 3D mesh (set ground and building heights)"
  {
    Timer timer("Step 3.5: Smooth 3D mesh (set ground and building heights)");
    LaplacianSmoother::SmoothMesh3D(mesh, cityModel, dtm, topHeight, true);
    Info(mesh);
  }

  // Write data for debugging and visualization
  if (p.WriteVTK)
  {
    Surface3D boundary;
    MeshProcessor::ExtractBoundary3D(boundary, mesh);
    VTK::Write(mesh, dataDirectory + "Step35Mesh.vtu");
    VTK::Write(boundary, dataDirectory + "Step35Boundary.vtu");
  }

  // Extract boundary of final mesh
  Surface3D boundary;
  MeshProcessor::ExtractBoundary3D(boundary, mesh);

  // Extract surface excluding top and sides
  Surface3D surface;
  MeshProcessor::ExtractOpenSurface3D(surface, boundary);

  // Write to file
  if (p.WriteJSON)
  {
    JSON::Write(mesh, dataDirectory + "CityMesh.json", origin);
    JSON::Write(surface, dataDirectory + "CitySurface.json", origin);
  }

  // Write data for debugging and visualization
  if (p.WriteVTK)
  {
    VTK::Write(mesh, dataDirectory + "CityMesh.vtu");
    VTK::Write(surface, dataDirectory + "CitySurface.vtu");
  }
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
  Parameters p;
  JSON::Read(p, argv[1]);
  Info(p);

  // Set data directory
  const std::string dataDirectory{p.DataDirectory + "/"};

  // Read city model
  CityModel cityModel;
  JSON::Read(cityModel, dataDirectory + "CityModel.json");
  Info(cityModel);

  // Read elevation model (only DTM is used)
  GridField2D dtm;
  JSON::Read(dtm, dataDirectory + "DTM.json");
  Info(dtm);

  // Generate surface meshes (non-matching, used for visualization)
  if (p.GenerateSurfaceMeshes)
  {
    GenerateSurfaceMeshes(cityModel, dtm, p);
  }

  // Generate volume meshes (matching, used for simulation)
  if (p.GenerateVolumeMeshes)
  {
    GenerateVolumeMeshes(cityModel, dtm, p);
  }

  // Report timings
  Timer::Report("dtcc-generate-mesh");

  return 0;
}
