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
  size_t xsize, ysize, bands = 0;
  BoundingBox2D bounds{};
  std::vector<GridField> values{};

  GeoRaster() {}

  double operator()(const Point2D& p, size_t band = 1) const {
    if (band > bands)
    {
      throw std::runtime_error("Raster only has " + str(bands) + " bands");
    }
    return values[band - 1].nearest(p);
  }

  double interpolate(const Point2D &p, size_t band = 1)
  {
    if (band > bands)
    {
      throw std::runtime_error("Raster only has " + str(bands) + " bands");
    }
    return values[band - 1](p);
  }
};
} // namespace DTCC_BUILDER

#endif
