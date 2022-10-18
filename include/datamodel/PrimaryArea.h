// Copyright (C) 2020 Anton J Olsson
// Licensed under the MIT License

#ifndef CORE_PRIMARYAREA_H
#define CORE_PRIMARYAREA_H

#include "BaseArea.h"
#include <Logging.h>
#include <vector>
namespace DTCC
{

/// Representation of a primary area.
class PrimaryArea : public Printable
{
public:
  /// Parent ID
  std::string DistrictAreaID;
  /// Primary area's name
  std::string Name;
  /// Base area constituting the primary area
  std::vector<BaseArea> BaseAreas;
  /// Base area's footprint
  Polygon Footprint;
  /// Base area's ID
  std::string AreaID;

  /// Pretty-print base area.
  /// \return Pretty-print string
  std::string __str__() const override
  {
    std::string baseAreaIDs;
    for (size_t i = 0; i < BaseAreas.size(); ++i)
      baseAreaIDs += (i == 0 ? "" : ", ") + BaseAreas[i].AreaID;
    return "Primary area with name " + Name + ", area ID " + AreaID +
           ", district area ID " + DistrictAreaID + " and base area(s) " +
           baseAreaIDs;
  }
};

} // namespace DTCC

#endif // CORE_PRIMARYAREA_H
