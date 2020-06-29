// Copyright (C) 2020 Dag Wästberg
// Licensed under the MIT License

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <gdal/cpl_conv.h>
#include <gdal/gdal_priv.h>

#include "Color.h"
#include "CSV.h"
#include "GeoRaster.h"
#include "RASTER.h"
#include "LAS.h"
#include "Logging.h"
#include "Point.h"
#include "PointCloud.h"
#include "PointCloudProcessor.h"

using namespace DTCC;

int main(int argc, char *argv[])
{

  PointCloud pointCloud;

  const std::string lasFile =
      "/home/dtcc/core/data/colorize_points/09B003_639_35_5025.las";
  const std::string tifFile =
      "/home/dtcc/core/data/colorize_points/1200_639_18_50_2018.tif";
  std::vector<int> cls = {1};
  BoundingBox2D bbox = BoundingBox2D(Point2D(182200,6395200),Point2D(182500,6395500));
  LAS::Read(pointCloud, lasFile, cls,bbox);

  std::cout << str(pointCloud.BoundingBox) << std::endl;
  GeoRaster orto;
  RASTER::Read(orto, tifFile,{1,2,3});
  
  std::cout << "XSize " << orto.XSize << " YSize " << orto.YSize << std::endl;
 
  PointCloudProcessor::ColorFromImage(pointCloud,orto);

  size_t i = 0;
  for (auto c: pointCloud.Colors) {
    std::cout << str(c) << std::endl;
    if (i == 10)
      break;
    i++;
  }

  LAS::Write(pointCloud,"/home/dtcc/core/data/colorize_points/testPC.las");

}