//
// Created by Anton Olsson on 2020-11-18.
//

#ifndef CORE_PRIMARYAREA_H
#define CORE_PRIMARYAREA_H

#include "BaseArea.h"
#include <Logging.h>
#include <vector>
namespace DTCC
{
class PrimaryArea : public Printable
{
public:
  std::string DistrictAreaID;
  std::string Name;
  std::vector<BaseArea> BaseAreas;
  Polygon Footprint;
  std::string AreaID;

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
