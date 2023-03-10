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
#include "MeshIO.h"
#include "MeshProcessor.h"
#include "OBJ.h"
#include "ParameterProcessor.h"
#include "Parameters.h"
#include "Utils.h"
#include "VTK.h"
#include "VertexSmoother.h"

#include "LAS.h"
#include "SHP.h"

using namespace DTCC_BUILDER;

void Help() { error("Usage: dtcc-generate-mesh Parameters.json"); }

// Generate surface meshes (non-matching, used for visualization)
void GenerateSurfaceMeshes(const CityModel &cityModel,
                           const GridField2D &dtm,
                           const Parameters &p)
{
  // Get data directory
  const std::string dataDirectory = p["DataDirectory"];

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

  if (p["WriteSTL"])
  {
    MeshIO::Write(groundSurface, dataDirectory + "GroundSurface.stl", "stl",
                  false);
    MeshIO::Write(buildingSurface, dataDirectory + "BuildingSurface.stl", "stl",
                  false);
  }

  if (p["WriteOBJ"])
  {
    MeshIO::Write(groundSurface, dataDirectory + "GroundSurface.obj", "obj",
                  true);
    MeshIO::Write(buildingSurface, dataDirectory + "BuildingSurface.obj", "obj",
                  true);
  }
}

// Generate volume meshes (matching, used for simulation)
void GenerateVolumeMeshes(CityModel &cityModel,
                          const GridField2D &dtm,
                          const Parameters &p)
{
  // Get directories
  const std::string dataDirectory = p["DataDirectory"];
  const std::string outputDirectory = p["OutputDirectory"];

  // Get origin (for serialization purposes)
  Point2D origin({p["X0"], p["Y0"]});

  // Step 1: Generate city model (and elevation model).
  // This step is handled by dtcc-generate-citymodel and
  // we assume that the data has already been generated.

  // Step 2.1: Merge building footprints
  {
    Timer timer("Step 2.1: Merge building footprints");
    CityModelGenerator::SimplifyCityModel(cityModel, dtm.Grid.BoundingBox,
                                          p["MinBuildingDistance"],
                                          p["MinVertexDistance"]);
    info(cityModel);
  }

  // Step 2.2: Clean building footprints
  {
    Timer timer("Step 2.2: Clean building footprints");
    CityModelGenerator::CleanCityModel(cityModel, p["MinVertexDistance"]);
    info(cityModel);
  }

  // Step 2.3: Compute building heights
  {
    Timer timer("Step 2.3: Compute building heights");
    CityModelGenerator::ComputeBuildingHeights(
        cityModel, dtm, p["GroundPercentile"], p["RoofPercentile"]);
    info(cityModel);
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
    info(mesh2D);
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
    info(mesh);
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
    LaplacianSmoother::SmoothMesh3D(mesh, cityModel, dtm, topHeight, false,
                                    false);
    info(mesh);
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
    info(mesh);
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
    LaplacianSmoother::SmoothMesh3D(mesh, cityModel, dtm, topHeight, true,
                                    p["WriteMatrix"]);
    info(mesh);
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
  if (p["WriteJSON"] || p["WriteVTK"] || p["WriteSTL"] || p["WriteOBJ"])
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

  if (p["WriteSTL"])
  {
    MeshIO::Write(surface, outputDirectory + "CitySurface.stl", "stl");
  }

  if (p["WriteOBJ"])
  {
    MeshIO::Write(surface, outputDirectory + "CitySurface.obj", "obj", true);
  }
}

int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (CommandLine::HasOption("-h", argc, argv))
  {
    Help();
    return 0;
  }
  Parameters p = ParameterProcessor::ProcessArgs(argc, argv);

  // Get directories
  const std::string dataDirectory = p["DataDirectory"];
  const std::string outputDirectory = p["OutputDirectory"];
  info("Loding data from directory: " + dataDirectory);
  info("Saving data to directory:   " + outputDirectory);

  // Read city model
  CityModel cityModel;
  JSON::Read(cityModel, outputDirectory + "CityModel.json");
  info(cityModel);

  // Read elevation model (only DTM is used)
  GridField2D dtm;
  JSON::Read(dtm, outputDirectory + "DTM.json");
  info(dtm);

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
  const std::string prefix = outputDirectory + "/dtcc-generate-mesh";
  Timer::Report("Timings for dtcc-generate-mesh", prefix + "-timings.json");
  JSON::Write(p, prefix + "-parameters.json", 4);
  info(p);

  return 0;
}
