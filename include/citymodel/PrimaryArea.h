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
class PrimaryArea
{
public:
  size_t DistrictAreaID;
  std::string Name;
  std::vector<BaseArea> BaseAreas;
  Polygon Footprint;
  size_t AreaID;
};
} // namespace DTCC

#endif // CORE_PRIMARYAREA_H
