// Copyright (C) 2021 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_CITY_MODEL_PROCESSOR_H
#define DTCC_CITY_MODEL_PROCESSOR_H

#include <vector>

#include "Geometry.h"

#include "datamodel/Building.h"
#include "datamodel/CityModel.h"

namespace DTCC
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
    info("CityModelProcessing: Filtered " + str(numFiltered) + " buildings");
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
    info("CityModelProcessing: Filtered " + str(numFiltered) + " buildings");
  }
};

} // namespace DTCC
#endif
