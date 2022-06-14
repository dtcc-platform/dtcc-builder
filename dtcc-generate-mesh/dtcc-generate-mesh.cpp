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
#include "ParameterProcessor.h"
#include "Parameters.h"
#include "Utils.h"
#include "VTK.h"
#include "VertexSmoother.h"

using namespace DTCC;

void Help() { Error("Usage: dtcc-generate-simulation-mesh Parameters.json"); }

// Generate surface meshes (non-matching, used for visualization)
void GenerateSurfaceMeshes(const CityModel &cityModel,
                           const GridField2D &dtm,
                           const Parameters &p)
{
  // Get data directory

  std::string dataDirectory = (std::string)p["DataDirectory"];
  std::string outputDirectory = (std::string)p["OutputDirectory"];
  // Get origin (for serialization purposes)
  Point2D origin({p["X0"], p["Y0"]});

  // Generate mesh for ground and buildings
  Surface3D groundSurface;
  std::vector<Surface3D> buildingSurfaces;
  MeshGenerator::GenerateSurfaces3D(groundSurface, buildingSurfaces, cityModel,
                                    dtm, p["MeshResolution"]);

  // Merge building surfaces
  Surface3D buildingSurface;
  MeshProcessor::MergeSurfaces3D(buildingSurface, buildingSurfaces);

  // Write to file
  if (p["WriteJSON"])
  {
    JSON::Write(groundSurface, dataDirectory + "GroundSurface.json", origin);
    JSON::Write(buildingSurface, dataDirectory + "BuildingSurface.json",
                origin);
  }

  // Write data for debugging and visualization
  if (p["WriteVTK"])
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
  // Get data directory
  std::string dataDirectory = (std::string)p["DataDirectory"];
  std::string outputDirectory = (std::string)p["OutputDirectory"];
  // Get origin (for serialization purposes)
  Point2D origin({p["X0"], p["Y0"]});

  // Step 1: Generate city model (and elevation model).
  // This step is handled by dtcc-generate-citymodel and
  // we assume that the data has already been generated.

  // Step 2.1: Merge building footprints
  {
    Timer timer("Step 2.1: Merge building footprints");
    CityModelGenerator::SimplifyCityModel(cityModel, p["MinBuildingDistance"],
                                          p["MinVertexDistance"]);
    Info(cityModel);
  }

  // Step 2.2: Clean building footprints
  {
    Timer timer("Step 2.2: Clean building footprints");
    CityModelGenerator::CleanCityModel(cityModel, p["MinVertexDistance"]);
    Info(cityModel);
  }

  // Step 2.3: Compute building heights
  {
    Timer timer("Step 2.3: Compute building heights");
    CityModelGenerator::ComputeBuildingHeights(
        cityModel, dtm, p["GroundPercentile"], p["RoofPercentile"]);
    Info(cityModel);
  }

  // Write JSON
  if (p["WriteJSON"])
  {
    JSON::Write(cityModel, dataDirectory + "CityModelSimple.json", origin);
  }

  // Step 3.1: Generate 2D mesh
  Mesh2D mesh2D;
  {
    Timer timer("Step 3.1: Generate 2D mesh");
    MeshGenerator::GenerateMesh2D(mesh2D, cityModel, dtm.Grid.BoundingBox,
                                  p["MeshResolution"]);
    Info(mesh2D);
  }

  // Write VTK
  if (p["WriteVTK"])
  {
    VTK::Write(mesh2D, dataDirectory + "Step31Mesh.vtu");
  }

  // Step 3.2: Generate 3D mesh (layer 3D mesh)
  size_t numLayers{};
  Mesh3D mesh;
  {
    Timer timer("Step 3.2: Generate 3D mesh");
    numLayers = MeshGenerator::GenerateMesh3D(mesh, mesh2D, p["DomainHeight"],
                                              p["MeshResolution"]);
    Info(mesh);
  }

  // Write VTK
  if (p["WriteVTK"])
  {
    Surface3D boundary;
    MeshProcessor::ExtractBoundary3D(boundary, mesh);
    VTK::Write(mesh, outputDirectory + "Step32Mesh.vtu");
    VTK::Write(boundary, dataDirectory + "Step32Boundary.vtu");
  };

  // Step 3.3: Smooth 3D mesh (set ground height)
  double topHeight{};
  {
    Timer timer("Step 3.3: Smooth 3D mesh");
    topHeight = dtm.Mean() + static_cast<double>(p["DomainHeight"]);
    LaplacianSmoother::SmoothMesh3D(mesh, cityModel, dtm, topHeight, false);
    Info(mesh);
  }

  // Write VTK
  if (p["WriteVTK"])
  {
    Surface3D boundary;
    MeshProcessor::ExtractBoundary3D(boundary, mesh);
    VTK::Write(mesh, outputDirectory + "Step33Mesh.vtu");
    VTK::Write(boundary, outputDirectory + "Step33Boundary.vtu");
  }

  // Step 3.4: Trim 3D mesh (remove building interiors)
  {
    Timer timer("Step 3.4: Trim 3D mesh");
    MeshGenerator::TrimMesh3D(mesh, mesh2D, cityModel, numLayers);
    Info(mesh);
  }

  // Write VTK
  if (p["WriteVTK"])
  {
    Surface3D boundary;
    MeshProcessor::ExtractBoundary3D(boundary, mesh);
    VTK::Write(mesh, outputDirectory + "Step34Mesh.vtu");
    VTK::Write(boundary, outputDirectory + "Step34Boundary.vtu");
  }

  // Step 3.5: Smooth 3D mesh (set ground and building heights)"
  {
    Timer timer("Step 3.5: Smooth 3D mesh");
    LaplacianSmoother::SmoothMesh3D(mesh, cityModel, dtm, topHeight, true);
    Info(mesh);
  }

  // Write VTK
  if (p["WriteVTK"])
  {
    Surface3D boundary;
    MeshProcessor::ExtractBoundary3D(boundary, mesh);
    VTK::Write(mesh, outputDirectory + "Step35Mesh.vtu");
    VTK::Write(boundary, outputDirectory + "Step35Boundary.vtu");
  }

  // Extract boundary of final mesh
  Surface3D boundary;
  Surface3D surface;
  if (p["WriteJSON"] || p["WriteVTK"])
  {
    MeshProcessor::ExtractBoundary3D(boundary, mesh);
    MeshProcessor::ExtractOpenSurface3D(surface, boundary);
  }

  // Write JSON
  if (p["WriteJSON"])
  {
    JSON::Write(mesh, outputDirectory + "CityMesh.json", origin);
    JSON::Write(surface, outputDirectory + "CitySurface.json", origin);
  }

  // Write VTK
  if (p["WriteVTK"])
  {
    VTK::Write(mesh, outputDirectory + "CityMesh.vtu");
    VTK::Write(surface, outputDirectory + "CitySurface.vtu");
  }
}

int main(int argc, char *argv[])
{
  // Check command-line arguments
  std::string dataDirectory;
  std::string outputDirectory;

  if (CommandLine::HasOption("-h", argc, argv))
  {
    Help();
    return 0;
  }
  Parameters p = ParameterProcessor::ProcessArgs(argc, argv);
  dataDirectory = (std::string)p["DataDirectory"];
  outputDirectory = (std::string)p["OutputDirectory"];

  Info("Loding from Data directory: " + dataDirectory);
  Info("Saving to Output directory: " + outputDirectory);

  // Read city model
  CityModel cityModel;
  JSON::Read(cityModel, outputDirectory + "CityModel.json");
  Info(cityModel);

  // Read elevation model (only DTM is used)
  GridField2D dtm;
  JSON::Read(dtm, outputDirectory + "DTM.json");
  Info(dtm);

  // Generate surface meshes (non-matching, used for visualization)
  if (p["GenerateSurfaceMeshes"])
  {
    GenerateSurfaceMeshes(cityModel, dtm, p);
  }

  // Generate volume meshes (matching, used for simulation)
  if (p["GenerateVolumeMeshes"])
  {
    GenerateVolumeMeshes(cityModel, dtm, p);
  }

  // Report timings and parameters
  Timer::Report("dtcc-generate-mesh", outputDirectory);
  Info(p);

  return 0;
}
