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
};

} // namespace DTCC

#endif // CORE_DISTRICT_H
