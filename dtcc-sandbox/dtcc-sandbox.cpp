// Copyright (C) 2020-2021 Anders Logg, Anton J Olsson
// Licensed under the MIT License

#include <string>
#include <vector>

// Needs to come before JSON (nlohmann) include because of sloppy
// namespacing in VTK (typedef detail)...
#include "VTK.h"

// DTCC includes
#include "CityModelGenerator.h"
#include "CommandLine.h"
#include "ElevationModelGenerator.h"
#include "Geometry.h"
#include "Grid.h"
#include "GridField.h"
#include "JSON.h"
#include "LAS.h"
#include "Logging.h"
#include "Parameters.h"
#include "Polygon.h"
#include "SHP.h"
#include "Timer.h"
#include "VertexSmoother.h"

using namespace DTCC;

void Help() { Error("Usage: dtcc-generate-citymodel Parameters.json"); }

void printSet(const std::set<size_t>& setToPrint)
{
  for(auto It = setToPrint.cbegin(); It != setToPrint.cend(); It++)
  {
    std::cout<<*It<<" ";
  }
  std::cout<<std::endl;
}

std::vector<std::vector<size_t>>
getNeighbors(size_t iId,
             const Grid2D &iGrid,
             const std::vector<Point2D> &iCenter,
             const std::vector<std::vector<size_t>> &iBins)
{
  std::vector<std::vector<size_t>> oNeighbors;
  size_t bin = iGrid.Point2Index(iCenter[iId]);
  std::vector<size_t> indices;
  oNeighbors.push_back(iBins[bin]);
  iGrid.Index2Boundary8(bin, indices);
  // e.g. 3,4,5
  for (size_t i = 0; i < indices.size(); i++)
  {
    // std::cout<<indices[i]<<std::endl;
    oNeighbors.push_back(iBins[indices[i]]);
  }
  return oNeighbors;
}

std::vector<std::set<size_t>>
getNeighborsSet(size_t iId,
                const Grid2D &iGrid,
                const std::vector<Point2D> &iCenter,
                const std::vector<std::set<size_t>> &iBins)
{

  std::vector<std::set<size_t>> oNeighbors;
  size_t bin = iGrid.Point2Index(iCenter[iId]);
  std::vector<size_t> indices;
  std::set<size_t> set;
  oNeighbors.push_back(iBins[bin]);
  iGrid.Index2Boundary8(bin, indices);
  // e.g. 3,4,5
  for (size_t i = 0; i < indices.size(); i++)
  {
    // set.insert(iBins[indices[i]]);

    oNeighbors.push_back(iBins[indices[i]]);    

    // set.clear();
    // std::cout<<indices[i]<<std::endl;
    // oNeighbors.push_back(iBins[indices[i]]);
  }
  return oNeighbors;
}

// std::set<size_t>
// getNeighborsSet(size_t id, const Grid2D& inGrid, const std::vector<Point2D>& inCenter, const std::vector<std::set<size_t>> bins, double distance)
// {
//   std::set<size_t> neighbors;
//   size_t bin = inGrid.Point2Index(inCenter[id]);
//   std::set<size_t> set;
//   std::vector<size_t> indices;
//   inGrid.Index2Boundary(bin,indices);
//   std::cout<<"Building:" << id << " is on the same bin with:"<<std::endl;
//   for(size_t i=0;i<indices.size();i++)
//   {
//     std::cout<<indices[i]<<std::endl;
//     neighbors.insert(bins[indices[i]]);
//   }
  
//   return neighbors;
// }

std::set<size_t>
getNeighborsSet(size_t footprintId, const std::vector<Polygon>& footprints, double distance)
{
  std::set<size_t> neighbors;
  const Polygon pivot = footprints[footprintId];
  for(size_t i=0;i<footprints.size();i++)
  {
    //Skip ourselves
    if(footprintId==i) continue;
    double dist2 = Geometry::SquaredDistance2D(pivot,footprints[i]);
    if(dist2<distance)
    {
      //std::cout<<"building: "<< footprintId<<" is neighbors with:" << i << " with a distance of:"<<dist2<<std::endl;
      neighbors.insert(i);
    }

  }

  std::cout<<"Building: " << footprintId<<" is neighbors with:"<<std::endl;
  printSet(neighbors);

  return neighbors;
}



bool canMergeFootprints(const Polygon& footprintToCheck, const std::vector<Polygon>& footprints, const std::set<size_t>& neighbors, size_t& mergeIndex, double tol)
{
  //const Polygon& a = footprints[startFootprintIndex];
  for (auto it = neighbors.cbegin(); it != neighbors.cend(); it++)
  {
    const Polygon& temp = footprints[*it];
    //const double d2 = Geometry::SquaredDistance2D(a,temp);
    const double d2 = Geometry::SquaredDistance2D(footprintToCheck, temp);
    if(d2<tol)
    {
      //std::cout<<"Can merge footprints: " << footprintToCheck<< " with: "<< *it << std::endl;
      std::cout<<"Can merge footprint:"<<*it<< " with a distance of:"<<d2<<std::endl;
      mergeIndex=*it;
      return true;
    }
  }
  std::cout<<"Can't merge any more footprints!"<<std::endl;
  return false;
}

void mergeFootprints(const size_t& startFootprintIndex,const std::vector<Polygon>& footprints, std::set<size_t> neighbors, double tol)
{
  size_t mergeIndex=0;
  Polygon pivotFootprint = footprints[startFootprintIndex];
  while(canMergeFootprints(pivotFootprint, footprints, neighbors, mergeIndex, tol))
  {
    auto neighborToDelete = neighbors.find(mergeIndex);
    if(mergeIndex == startFootprintIndex)
    {
      neighbors.erase(neighborToDelete);
      std::cout<<"Tried to merge footprint with itself, skipping..."<<std::endl;
      continue;
    }

    Polygon footprintB = footprints[*neighborToDelete];
    
    pivotFootprint = Polyfix::MergePolygons(pivotFootprint, footprintB, tol);

    neighbors.erase(neighborToDelete);

    std::cout<<"Erased footprint:"<<*neighborToDelete<<std::endl;
    std::cout<<"Printin remaining set..."<<std::endl;
    printSet(neighbors);
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
  const std::string modelName = Utils::GetFilename(argv[1], true);

  // Get data directory
  std::string dataDirectory = p["DataDirectory"];
  dataDirectory += "/";

  // Start timer
  Timer timer("Step 1: Generate city model");

  // Read property map
  std::vector<Polygon> footprints;
  std::vector<std::string> UUIDs;
  std::vector<int> entityIDs;
  SHP::Read(footprints, dataDirectory + "PropertyMap.shp", &UUIDs, &entityIDs);
  Info("Loaded " + str(footprints.size()) + " building footprints");

  // Set bounding box
  BoundingBox2D bbox;
  Point2D O;
  if (p["AutoDomain"])
  {
    bbox = BoundingBox2D(footprints, p["DomainMargin"]);
    Info("Bounding box of footprints: " + str(bbox));
    BoundingBox2D lasBBox;
    LAS::BoundsDirectory(lasBBox, dataDirectory);
    Info("Bounding box of point cloud: " + str(lasBBox));
    bbox.Intersect(lasBBox);
    O = bbox.P;
    p["X0"] = O.x;
    p["Y0"] = O.y;
    p["XMin"] = 0.0;
    p["YMin"] = 0.0;
    p["XMax"] = bbox.Q.x - bbox.P.x;
    p["YMax"] = bbox.Q.y - bbox.Q.y;
  }
  else
  {
    O = Point2D(p["X0"], p["Y0"]);
    const double xMin = p["XMin"];
    const double xMax = p["XMax"];
    const double yMin = p["YMin"];
    const double yMax = p["YMax"];
    const Point2D P{O.x + xMin, O.y + yMin};
    const Point2D Q{O.x + xMax, O.y + yMax};
    bbox = BoundingBox2D(P, Q);
  }

  // Check size of bounding box
  Info("Bounding box of city model: " + str(bbox));
  if (bbox.Area() < 100.0)
  {
    Error("Domain too small to generate a city model");
    return 1;
  }
  std::vector<double> diameters;
  std::vector<Point2D> centers;
  double maxdiam = 0;
  for (int i = 0; i < footprints.size(); i++)
  {
    centers.push_back(Geometry::PolygonCenter2D(footprints[i]));
    diameters.push_back(Geometry::PolygonRadius2D(footprints[i], centers[i]));
    std::cout << str(centers[i]) << " " << str(diameters[i]) << " " << i
              << std::endl;
    if (diameters[i] > maxdiam)
    {
      maxdiam = diameters[i];
    }
  }
  std::cout << maxdiam << std::endl;
  double xlength = bbox.Q.x - bbox.P.x;
  double ylength = bbox.Q.y - bbox.P.y;

  size_t x = (bbox.Q.x - bbox.P.x) / (2 * maxdiam);
  size_t y = (bbox.Q.y - bbox.P.y) / (2 * maxdiam);

  Info("X Length:" + str(xlength));
  Info("Y Length:" + str(ylength));

  Info("X:" + str(x));
  Info("Y:" + str(y));

  size_t temp;
  std::cin>>temp;
  Grid2D grid(bbox, x, y);
  std::vector<std::set<size_t>> neighbors; //,setbins;
  std::vector<std::vector<size_t>> bins;
  std::vector<std::set<size_t>> setbins;
  bins.reserve(x * y);
  setbins.reserve(x * y);
  bins = {std::vector<size_t>()};
  for (size_t i = 0; i < x * y; i++)
  {
    bins.push_back(std::vector<size_t>());
    setbins.push_back(std::set<size_t>());
  }
  // std::set<size_t> set;
  for (int i = 0; i < footprints.size(); i++)
  {
    Info("Building " + str(i) + " " + str(grid.Point2Index(centers[i])));
    // ins[grid.Point2Index(centers[i])]=i;
  }

  // std::vector<size_t> indices;

  /*   grid.Index2Boundary8(15, indices);
    for (int i=0;i<indices.size();i++)
    {
    set.insert(i);
    neighbors.push_back(set);
    std::cout<<indices[i]<<std::endl;
    }
    std::string dummy;
    std::getline(std::cin, dummy);
    indices.clear();
     grid.Index2Boundary(15, indices);
    for (int i=0;i<indices.size();i++)
    {
    set.insert(i);
    neighbors.push_back(set);
    std::cout<<indices[i]<<std::endl;
    }
    std::getline(std::cin, dummy); */
  std::set<size_t> myset;
  Info("I am here 209");
  for (size_t i = 0; i < footprints.size(); i++)
  {
    Info("Building " + str(i) + " " + str(grid.Point2Index(centers[i])));
    size_t vertexno = grid.Point2Index(centers[i]);
    // bins[grid.Point2Index(centers[i])].push_back(i);
    // myset.insert(i);
    std::cout << "Vertexno " << vertexno << std::endl;
    setbins[vertexno].insert(i);
    // ins[grid.Point2Index(centers[i])]=i;
  }

  for (size_t i = 0; i < footprints.size(); i++)
  {
    Info("Building " + str(i) + " " + str(grid.Point2Index(centers[i])));
    size_t vertexno = grid.Point2Index(centers[i]);
    // bins[grid.Point2Index(centers[i])].push_back(i);
    bins[vertexno].push_back(i);
    // ins[grid.Point2Index(centers[i])]=i;
  }
  Info("Bins size is " + str(bins.size()));
  size_t sum = 0;
  for (size_t i = 0; i < bins.size(); i++)
  {
    sum = sum + bins[i].size();
    std::cout << "Bin " << i << " has size " << bins[i].size() << " and total "
              << sum << std::endl;
  }

  for (size_t i = 0; i < bins.size(); i++)
  {
    Info("Bin " + str(i) + " has size " + str(bins[i].size()) +
         " and buildings: ");
    for (size_t j = 0; j < bins[i].size(); j++)
    {
      Info(str(bins[i][j]) + " " + str(j));
    }
  }

  /* for (size_t i = 0; i < bins.size(); i++)
  {
    Info("Near cells of bin " + str(i) + " are:");
    for (size_t j = 0; j < bins[i].size(); j++)
    {
      std::vector<size_t> indices;
      grid.Index2Boundary8(bins[i][j], indices);
      for (size_t k = 0; k < indices.size(); k++)
        Info(str(indices[k]) + " " + str(i) + " " + str(j));
    }
    indices.clear();
  } */

  // std::cout << "I am here 261" << std::endl;
  // std::vector<std::vector<size_t>> nearest =
  //     getNeighbors(0, grid, centers, bins);
  // std::vector<std::set<size_t>> nearestset =
  //     getNeighborsSet(0, grid, centers, setbins);

  // std::cout << nearest.size() << std::endl;
  // for (size_t i = 0; i < nearest.size(); i++)
  // {
  //   for (size_t j = 0; j < nearest[i].size(); j++)
  //   {
  //     std::cout << "Nearest buildings " << nearest[i][j] << std::endl;
  //   }
  // }

  //neighbor print....
  std::cout<<"Printing neighbors of buildings with getNeighborsSet...."<<std::endl;
  for (size_t i = 0; i < footprints.size(); i++)
  {
    std::vector<std::set<size_t>> neighborsVecSet =
        getNeighborsSet(i, grid, centers, setbins);
    std::cout<<"-------------"<<std::endl;
    Info("Neighbors of building " + str(i) + " are ");
    //Info(footprints[i]);
    for (size_t j = 0; j < neighborsVecSet.size(); j++)
    {
      for (auto it = neighborsVecSet[j].cbegin(); it != neighborsVecSet[j].cend(); it++)
      {
        std::cout << *it << ' ';
      }
      std::cout << std::endl;
    }
    // for (size_t j = 0; j < neighborsVecSet.size(); j++)
    // {
    //   std::cout<<"Trying to merge set #"<<j<<std::endl;
    //   mergeFootprints(i,footprints,neighborsVecSet[j],p["MinBuildingDistance"]);
    // }
  }

  std::cout<<"Printing neighbors of buildings with brute force...."<<std::endl;
  std::cout<<"bb X:"<<grid.XStep<<" , "<<grid.YStep<<std::endl;
  size_t t;
  std::cin>>t;

  for(size_t i=0;i<footprints.size();i++)
  {
    //std::set<size_t> neighbors = getNeighborsSet(i,footprints,bbox.Area());
    //printSet(neighbors);

    //Brute force...
    //std::set<size_t> neighbors = getNeighborsSet(i,footprints,p["MinBuildingDistance"]);
    std::set<size_t> neighbors = getNeighborsSet(i,footprints,grid.XStep);

    //Bin usage
    //std::set<size_t> neighbors2 = getNeighborsSet(i, grid, centers, 250.0);
    std::cout<<"------"<<std::endl;
  }

  // double tolerance = p["MinBuildingDistance"];
  // std::queue<size_t> indices;
  // for(size_t i=0;i<footprints.size();i++)
  //   indices.push(i);
  
  // size_t NumMerged=0;
  // while(!indices.empty())
  // {
  //   //pop indice of next building to check
  //   const size_t i=indices.front();
  //   indices.pop();
  //   std::vector<std::set<size_t>> neighbors = getNeighborsSet(i, grid, centers, setbins);
  //   for(size_t j=0;j<neighbors.size();j++)
  //   {
  //     for(auto it = neighbors[j].cbegin();it!=neighbors[j].cend();it++)
  //     {
  //       // skip building itself
  //       if(i==j) continue;

  //       //if building[j] is empty skip

  //       //check distance for neighbors 
  //       const Polygon& Pi = footprints[i];
  //       const Polygon& Pj = footprints[neighbors[j]];
  //       const double d2 = Geometry::SquaredDistance2D(Pi,Pj);

  //       //Merge if distance is too small
  //       if(d2 < tol2)
  //       {
  //         //TODO: Call Merge Buildings here
  //         NumMerged++;
  //         indices.push(i);
  //       }
  //     }
  //   }
  // }

  /*   for (int i = 0; i < (x * y); i++)
    {
      Info("Bin " + str(i) + " has size " + str(bins[i].size()));
    }
    indices.clear();

    for (int i = 0; i < footprints.size(); i++)
    {
      size_t current = grid.Point2Index(centers[i]);
      Info("I am building " + str(i) + " I belong to bin " + str(current) +
           " and my neighbor cells are: ");
      grid.Index2Boundary8(current, indices);
      Info("And my neighbour buildings are:");

      for (int j = 0; j < indices.size(); j++)
      {
        std::cout << indices[j] << std::endl;
        for (int k = 0; k < bins[indices[j]].size(); k++)
        {
          std::cout << bins[indices[j]][k];
        }
      }

      indices.clear();
    }

    std::string dummy;
    std::getline(std::cin, dummy);
    size_t k = 0;
    for (int i = 0; i < x * y; i++)
      for (int j = 0; j < bins[i].size(); j++)
      {
        {
          Info("I am bin " + str(i) + " and I have buildings " +
  str(bins[i][j])); k++;
          // std::getline(std::cin, dummy);
          std::cout << "counter " << k << std::endl;
        }
      }
    std::getline(std::cin, dummy);
      for (int i = 0; i < footprints.size(); i++)
  {
      size_t cell=grid.Point2Index(centers[i]);
      Info("I am building "+str(i)+", I belong in bin " +str(cell));
      for (int j=0;j<bins[cell].size();j++)
      {
      Info("Buildings in same cell are buildings "+str(bins[cell][j]));
      }
      Info("Cells attached to mine are ");

      for (int j=0;j<bins[cell].size();j++)
      {
      Info("Buildings in same cell are buildings "+str(bins[cell][j]));
      }
  }

   */

  // Read point cloud (only points inside bounding box)
  /*
  PointCloud pointCloud;
  LAS::ReadDirectory(pointCloud, dataDirectory, bbox);

  // Check point cloud
  if (pointCloud.Empty())
    Error("Point cloud is empty. Check LiDaR quality or the X{0,Min,Max}, "
          "Y{0,Min,Max} values in Parameters.json");
  pointCloud.SetOrigin(O);
  Info(pointCloud);

  // Remove outliers from point cloud
  PointCloudProcessor::RemoveOutliers(pointCloud, p["OutlierMargin"]);

  // Generate DTM (excluding buildings and other objects)
  GridField2D dtm;
  ElevationModelGenerator::GenerateElevationModel(
      dtm, pointCloud, {2, 9}, p["ElevationModelResolution"]);
  Info(dtm);

  // Generate DSM (including buildings and other objects)
  GridField2D dsm;
  ElevationModelGenerator::GenerateElevationModel(
      dsm, pointCloud, {}, p["ElevationModelResolution"]);
  Info(dsm);

  // Smooth elevation model (only done for DTM)
  VertexSmoother::SmoothField(dtm, p["GroundSmoothing"]);

  // Generate raw city model
  CityModel cityModel;
  cityModel.Name = modelName;
  CityModelGenerator::GenerateCityModel(cityModel, footprints, UUIDs, entityIDs,
                                        bbox, p["MinBuildingDistance"]);
  cityModel.SetOrigin(O);
  Info(cityModel);

  // Clean city model and compute heights
  CityModelGenerator::CleanCityModel(cityModel, p["MinVertexDistance"]);
  CityModelGenerator::ExtractBuildingPoints(
      cityModel, pointCloud, p["GroundMargin"], p["OutlierMargin"]);
  CityModelGenerator::ComputeBuildingHeights(
      cityModel, dtm, p["GroundPercentile"], p["RoofPercentile"]);

  // Stop timer
  timer.Stop();

  // Write JSON
  if (p["WriteJSON"])
  {
    JSON::Write(dtm, dataDirectory + "DTM.json", O);
    JSON::Write(dsm, dataDirectory + "DSM.json", O);
    JSON::Write(cityModel, dataDirectory + "CityModel.json", O);
  }

  // Write VTK
  if (p["WriteVTK"])
  {
    VTK::Write(dtm, dataDirectory + "DTM.vts");
    VTK::Write(dsm, dataDirectory + "DSM.vts");
  }

  // Report timings and parameters
  Timer::Report("dtcc-generate-citymodel", dataDirectory);
  Info(p);
  */
  return 0;
}
