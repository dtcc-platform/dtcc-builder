// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#include <deque>
#include <dolfin.h>
#include <iostream>
#include <string>
#include <vector>

#include "CityModel.h"
#include "CommandLine.h"
#include "GridField.h"
#include "JSON.h"
#include "VTK.h"
#include "Mesh.h"
#include "MeshGenerator.h"
#include "Parameters.h"
#include "Logging.h"

using namespace DTCC;

void Help()
{
  Error("Usage: dtcc-generate-visualization-mesh Parameters.json");
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
  Parameters parameters;
  JSON::Read(parameters, argv[1]);
  Info(parameters);

  // Get data directory (add trailing slash just in case)
  const std::string dataDirectory = parameters.DataDirectory + "/";

  // Read city model data
  CityModel cityModel;
  JSON::Read(cityModel, dataDirectory + "CityModelClean.json");
  Info(cityModel);

  // Read terrain model
  GridField2D dtm{};
  JSON::Read(dtm, dataDirectory + "DTM.json");
  Info(dtm);

  // Generate 3D surfaces
  std::vector<Surface3D> surfaces = MeshGenerator::GenerateSurfaces3D(
      cityModel, dtm, dtm.Grid.BoundingBox.P.x, dtm.Grid.BoundingBox.P.y,
      dtm.Grid.BoundingBox.Q.x, dtm.Grid.BoundingBox.Q.y,
      parameters.MeshResolution, parameters.FlatGround);

  // Extract ground surface
  Surface3D groundSurface = surfaces[0];

  // FIXME: Consider moving this somewhere else as a utility for merging surfaces
  // Extract building surface as one common surface
  Progress("Merging building meshes...");
  Surface3D buildingSurface;
  size_t numPoints = 0;
  size_t numCells = 0;
  for (size_t i = 1; i < surfaces.size(); i++)
  {
    numPoints += surfaces[i].Vertices.size();
    numCells += surfaces[i].Faces.size();
  }
  buildingSurface.Vertices.resize(numPoints);
  buildingSurface.Faces.resize(numCells);
  size_t k = 0;
  size_t l = 0;
  for (size_t i = 1; i < surfaces.size(); i++)
  {
    for (size_t j = 0; j < surfaces[i].Faces.size(); j++)
    {
      Simplex2D c = surfaces[i].Faces[j];
      c.v0 += k;
      c.v1 += k;
      c.v2 += k;
      buildingSurface.Faces[l++] = c;
    }
    for (size_t j = 0; j < surfaces[i].Vertices.size(); j++)
      buildingSurface.Vertices[k++] = surfaces[i].Vertices[j];
  }

  // Write to files
  JSON::Write(groundSurface, dataDirectory + "VisualizationGroundMesh.json");
  VTK::Write(groundSurface, dataDirectory + "VisualizationGroundMesh.vtu");
  VTK::Write(buildingSurface, dataDirectory + "VisualizationBuildingMesh.vtu");

  // Report timings
  Timer::Report("dtcc-generate-visualization-mesh");

  return 0;
}
