// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <gdal/cpl_conv.h>
#include <gdal/gdal_priv.h>

#include "Color.h"
#include "GeoRaster.h"
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
      "/home/dtcc/core/data/colorize_points/09B003_639_35_5025_trees.las";
  const std::string tifFile =
      "/home/dtcc/core/data/colorize_points/1200_639_18_50_2018.tif";
  LAS::Read(pointCloud, lasFile);

  std::cout << str(pointCloud.BoundingBox) << std::endl;

  GeoRaster orto = GeoRaster(tifFile);
  
  std::cout << "XSize " << orto.XSize << " YSize " << orto.YSize << std::endl;
 
  PointCloudProcessor::ColorFromImage(pointCloud,orto);

  size_t i = 0;
  for (auto c: pointCloud.Colors) {
    std::cout << str(c) << std::endl;
    if (i == 10)
      break;
    i++;
  }

}