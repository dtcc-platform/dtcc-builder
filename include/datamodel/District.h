// Copyright (C) 2020 Anton J Olsson
// Licensed under the MIT License

#ifndef CORE_DISTRICT_H
#define CORE_DISTRICT_H

#include "PrimaryArea.h"
#include <Logging.h>
#include <Polygon.h>
namespace DTCC_BUILDER
{

/// Representation of a city district
class District : public Printable
{
public:
  /// District's footprint
  Polygon Footprint;
  /// District's name
  std::string Name;
  /// District's area ID
  std::string AreaID;
  /// Primary areas constituting the district
  std::vector<PrimaryArea> PrimaryAreas;

  /// Pretty-print district.
  /// \return Pretty-print string.
  std::string __str__() const override
  {
    std::string primAreaIDs;
    for (size_t i = 0; i < PrimaryAreas.size(); ++i)
    {
      primAreaIDs += (i == 0 ? "" : ", ") + PrimaryAreas[i].AreaID;
    }
    return "District with name " + Name + ", area ID " + AreaID +
           " and primary areas " + primAreaIDs;
  }
};

} // namespace DTCC_BUILDER

#endif // CORE_DISTRICT_H
