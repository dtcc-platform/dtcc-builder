// Copyright (C) 2021 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_CITY_MODEL_PROCESSOR_H
#define DTCC_CITY_MODEL_PROCESSOR_H

#include <vector>

#include "Geometry.h"

#include "Building.h"
#include "CityModel.h"

namespace DTCC_BUILDER
{

class CityModelProcessor
{
public:
  static void ErrorFilter(const CityModel &baseModel,
                          CityModel &outModel,
                          size_t errorFilter)
  {
    if (outModel.Name.empty())
    {
      outModel.Name = baseModel.Name;
    }
    size_t numFiltered = 0;
    for (auto &building : baseModel.Buildings)
    {
      if ((building.error & errorFilter) == 0)
      {
        outModel.Buildings.push_back(building);
      }
      else
      {
        numFiltered++;
      }
    }
    info("Filtered " + str(numFiltered) + " buildings");
  }

  static void BuildingFootprintFilter(const CityModel &baseModel,
                                      CityModel &outModel,
                                      double minArea)
  {
    if (outModel.Name.empty())
    {
      outModel.Name = baseModel.Name;
    }

    size_t numFiltered = 0;
    for (auto &building : baseModel.Buildings)
    {
      if (Geometry::PolygonArea(building.Footprint) >= minArea)
      {
        outModel.Buildings.push_back(building);
      }
      else
      {
        numFiltered++;
      }
    }
    info("Filtered " + str(numFiltered) + " buildings");
  }
};

} // namespace DTCC_BUILDER
#endif
