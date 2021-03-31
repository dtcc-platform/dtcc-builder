#include "JSON.h"
#include "Mesh.h"
#include <dolfin.h>
//#include "MeshField.h"

#include "FEniCS.h"

#include "CityModelGenerator.h"
#include "LaplacianSmoother.h"
#include "Logging.h"
#include "MeshGenerator.h"
#include "MeshProcessor.h"
#include "VTK.h"
//#include <vtkCellArray.h>
//#include <vtkCellData.h>
//#include <vtkFloatArray.h>
//#include <vtkPointData.h>
//#include <vtkPoints.h>
//#include <vtkSmartPointer.h>
//#include <vtkStructuredGrid.h>
//#include <vtkXMLStructuredGridWriter.h>
using namespace DTCC;
void Help() { Error("Usage: dtcc-generate-simulation-mesh Parameters.json"); }
void GenerateVolumeMeshes(CityModel &cityModel,
                          const GridField2D &dtm,
                          const Parameters &p)
{
  // Set data directory
  const std::string dataDirectory{p.DataDirectory + "/"};
  // Step 1: Generate city model (and elevation model).
  // This step is handled by dtcc-generate-citymodel and
  // we assume that the data has already been generated.
  // Step 2: Simplify city model (merge buildings)
  // CityModelGenerator::SimplifyCityModel(cityModel, p.MinBuildingDistance,
  //                                        p.MinVertexDistance);
  Info(cityModel);
  // Recompute building heights (seems better not to recompute heights).
  // If we recompute the heights, then we tend to get a height that is too low.
  // CityModelGenerator::ComputeBuildingHeights(cityModel, dtm,
  // p.GroundPercentile,
  //                                           p.RoofPercentile);
  // Step 3.1: Generate 2D mesh
  Mesh2D mesh2D;
  MeshGenerator::GenerateMesh2D(mesh2D, cityModel, dtm.Grid.BoundingBox,
                                p.MeshResolution);
  Info(mesh2D);
  // Write data for debugging and visualization
  if (p.Debug)
  {
    VTK::Write(mesh2D, dataDirectory + "Step31Mesh.vtu");
  }
  // Step 3.2: Generate 3D mesh (full)
  Mesh3D mesh;
  const size_t numLayers = MeshGenerator::GenerateMesh3D(
      mesh, mesh2D, p.DomainHeight, p.MeshResolution);
  Info(mesh);
  // Write data for debugging and visualization
  if (p.Debug)
  {
    Surface3D boundary;
    MeshProcessor::ExtractBoundary3D(boundary, mesh);
    VTK::Write(mesh, dataDirectory + "Step32Mesh.vtu");
    VTK::Write(boundary, dataDirectory + "Step32Boundary.vtu");
  };
  // Step 3.3: Smooth 3D mesh (apply DTM to ground)
  // const double topHeight = dtm.Mean() + p.DomainHeight;
  // LaplacianSmoother::SmoothMesh3D(mesh, cityModel, dtm, topHeight, false);
  Info(mesh);
  // Write data for debugging and visualization
  if (p.Debug)
  {
    Surface3D boundary;
    MeshProcessor::ExtractBoundary3D(boundary, mesh);
    VTK::Write(mesh, dataDirectory + "Step33Mesh.vtu");
    VTK::Write(boundary, dataDirectory + "Step33Boundary.vtu");
  }
  // Step 3.4: Trim 3D mesh (remove tets inside buildings)
  // MeshGenerator::TrimMesh3D(mesh, mesh2D, cityModel, numLayers);
  Info(mesh);
  // Write data for debugging and visualization
  if (p.Debug)
  {
    Surface3D boundary;
    MeshProcessor::ExtractBoundary3D(boundary, mesh);
    VTK::Write(mesh, dataDirectory + "Step34Mesh.vtu");
    VTK::Write(boundary, dataDirectory + "Step34Boundary.vtu");
  }
  // Step 3.5: Smooth 3D mesh (apply DTM to ground)
  // LaplacianSmoother::SmoothMesh3D(mesh, cityModel, dtm, topHeight, true);
  Info(mesh);
  // Write data for debugging and visualization
  if (p.Debug)
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
  // Get origin (for serialization purposes)
  Point2D origin({p.X0, p.Y0});
  // Write to file
  JSON::Write(mesh, dataDirectory + "CityMesh.json", origin);
  JSON::Write(surface, dataDirectory + "CitySurface.json", origin);
  // Write data for debugging and visualization
  if (p.Debug)
  {
    VTK::Write(mesh, dataDirectory + "CityMesh.vtu");
    VTK::Write(surface, dataDirectory + "CitySurface.vtu");
  }

  // convert to fenics mesh
  dolfin::Mesh d_mesh;
  FEniCS::ConvertMesh(mesh, d_mesh);
  // write/save to fenics mesh
  FEniCS::Write(d_mesh, dataDirectory + "test.xml");
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
  GridField2D dtm;
  // Point2D P={-150,-200};
  // Point2D Q={150,1250};
  // BoundingBox2D BBox={P,Q};
  // dtm.Grid.BoundingBox=BBox;
  const Point2D O{p.X0, p.Y0};
  const Point2D P{p.XMin + p.X0, p.YMin + p.Y0};
  const Point2D Q{p.XMax + p.X0, p.YMax + p.Y0};
  const BoundingBox2D bbox{P, Q};
  Info("Bounding box: " + str(bbox));
  CityModel cityModel;
  // const Point2D origin{0.0, 0.0};
  const auto resolution = p.ElevationModelResolution;
  dtm.Grid.XSize =
      (dtm.Grid.BoundingBox.Q.x - dtm.Grid.BoundingBox.P.x) / resolution + 1;
  dtm.Grid.YSize =
      (dtm.Grid.BoundingBox.Q.y - dtm.Grid.BoundingBox.P.y) / resolution + 1;
  dtm.Values.resize(dtm.Grid.XSize * dtm.Grid.YSize);
  std::fill(dtm.Values.begin(), dtm.Values.end(), 0.0);
  dtm.Grid.XStep = (dtm.Grid.BoundingBox.Q.x - dtm.Grid.BoundingBox.P.x) /
                   (dtm.Grid.XSize - 1);
  dtm.Grid.YStep = (dtm.Grid.BoundingBox.Q.y - dtm.Grid.BoundingBox.P.y) /
                   (dtm.Grid.YSize - 1);

  const double SCALE = 1; // millimeters
  const double D = 100.0 * SCALE;
  const double L = 50.00 * SCALE;
  const double H = 48.0 * SCALE;

  for (int i = 0; i < 30; i++)
  {
    // each even row has 15 houses and 14 for odd
    const int numHouses = (i % 2 == 0) ? 15 : 14;

    for (int j = 0; j < numHouses; j++)
    {

      // Make sure the houses are in a staggered formation
      const double shift = ((i % 2) + 1) * L - 15 * L;
      const Point2D position =
          Point2D(j * D + shift, (i + 0.5) * D - 2000.0 * SCALE);
      Building building =
          CityModelGenerator::GenerateBuilding(position, L, L, H, 0.0);
      cityModel.Buildings.push_back(building);
    }
  }
  JSON::Write(cityModel, dataDirectory + "FredrikCityModelRandom.json", O);
  JSON::Write(dtm, dataDirectory + "Dtm.json", O);

  Mesh2D mesh2D;
  MeshGenerator::GenerateMesh2D(mesh2D, cityModel, dtm.Grid.BoundingBox,
                                p.MeshResolution);
  Info(mesh2D);

  GenerateVolumeMeshes(cityModel, dtm, p);

  return EXIT_SUCCESS;
}
