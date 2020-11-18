//
// Created by Anton Olsson on 2020-11-18.
//

#ifndef CORE_BASEAREA_H
#define CORE_BASEAREA_H

#include "Property.h"
#include <vector>
namespace DTCC
{
class BaseArea
{
public:
  size_t AreaID;
  size_t PrimaryAreaID;
  std::vector<Property> Properties;
  std::vector<Building> Buildings;
};
} // namespace DTCC

#endif // CORE_BASEAREA_H
