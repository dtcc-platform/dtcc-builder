// Copyright (C) 2020 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_GEO_RASTER_H
#define DTCC_GEO_RASTER_H


#include "BoundingBox.h"
#include "GridField.h"
#include "Point.h"
#include "Logging.h"
namespace DTCC_BUILDER
{
class GeoRaster
{
public:
  size_t XSize, YSize, Bands = 0;
  BoundingBox2D Bounds{};
  std::vector<GridField> Values{};

  GeoRaster() {}

  double operator()(const Point2D& p, size_t band = 1) const {
    if (band > Bands) 
    {
      throw std::runtime_error("Raster only has " + str(Bands) + " bands");
    }
    return Values[band-1].Nearest(p);
  }

  double Interpolate(const Point2D& p, size_t band = 1) {
    if (band > Bands) 
    {
      throw std::runtime_error("Raster only has " + str(Bands) + " bands");
    }
    return Values[band-1](p);
  }
};
} // namespace DTCC_BUILDER

#endif
