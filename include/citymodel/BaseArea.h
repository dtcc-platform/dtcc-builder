//
// Created by Anton Olsson on 2020-11-18.
//

#ifndef CORE_BASEAREA_H
#define CORE_BASEAREA_H

#include "Property.h"
#include <vector>
namespace DTCC
{
class BaseArea : public Printable
{
public:
  size_t AreaID;
  size_t PrimaryAreaID;
  Polygon Footprint;
  std::vector<Property> Properties;
  std::vector<Building> Buildings;

  std::string __str__() const override
  {
    return "Base area with area ID " + std::to_string(AreaID) +
           ", primary area ID " + std::to_string(PrimaryAreaID) + ", " +
           std::to_string(Properties.size()) + " properties and " +
           std::to_string(Buildings.size()) + " buildings";
  }
};

} // namespace DTCC

#endif // CORE_BASEAREA_H
