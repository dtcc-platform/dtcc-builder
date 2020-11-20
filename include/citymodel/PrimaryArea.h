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
  size_t DistrictAreaID;
  std::string Name;
  std::vector<BaseArea> BaseAreas;
  Polygon Footprint;
  size_t AreaID;

  std::string __str__() const override
  {
    std::string baseAreaIDs;
    for (size_t i = 0; i < BaseAreas.size(); ++i)
      baseAreaIDs += (i == 0 ? "" : ", ") + str(BaseAreas[i].AreaID);
    return "Primary area with name " + Name + ", area ID " + str(AreaID) +
           ", district area ID " + std::to_string(DistrictAreaID) +
           " and base area(s) " + baseAreaIDs;
  }
};

} // namespace DTCC

#endif // CORE_PRIMARYAREA_H
