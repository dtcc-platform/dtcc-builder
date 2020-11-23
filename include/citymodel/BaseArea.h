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
  std::string AreaID;
  std::string PrimaryAreaID;
  Polygon Footprint;
  std::vector<Property> Properties;
  std::vector<Building> Buildings;

  std::string __str__() const override
  {
    return "Base area with area ID " + AreaID + ", primary area ID " +
           PrimaryAreaID + ", " + std::to_string(Properties.size()) +
           " properties and " + std::to_string(Buildings.size()) + " buildings";
  }
};

} // namespace DTCC

#endif // CORE_BASEAREA_H
