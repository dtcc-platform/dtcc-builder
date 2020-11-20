// Copyright (C) 2020 Anton J Olsson
// Licensed under the MIT License

#ifndef CORE_DISTRICT_H
#define CORE_DISTRICT_H

#include "PrimaryArea.h"
#include <Logging.h>
#include <Polygon.h>
namespace DTCC
{

class District : public Printable
{
public:
  Polygon Footprint;
  std::string Name;
  size_t AreaID;
  std::vector<PrimaryArea> PrimaryAreas;

  std::string __str__() const override
  {
    std::string primAreaIDs;
    for (size_t i = 0; i < PrimaryAreas.size(); ++i)
    {
      primAreaIDs +=
          (i == 0 ? "" : ", ") + std::to_string(PrimaryAreas[i].AreaID);
    }
    return "District with name " + Name + ", area ID " +
           std::to_string(AreaID) + " and primary areas " + primAreaIDs;
  }
};

} // namespace DTCC

#endif // CORE_DISTRICT_H
