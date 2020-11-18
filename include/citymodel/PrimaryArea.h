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
  size_t AreaID;
  size_t DistrictAreaID;
  std::string Name;
  std::vector<BaseArea> BaseAreas;
  std::vector<Polygon> Footprint;
};
} // namespace DTCC

#endif // CORE_PRIMARYAREA_H
