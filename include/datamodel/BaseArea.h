// Copyright (C) 2020 Anton J Olsson
// Licensed under the MIT License

#ifndef CORE_BASEAREA_H
#define CORE_BASEAREA_H

#include "Property.h"
#include <vector>
namespace DTCCBUILDER
{
/// Representation of a base area.
class BaseArea : public Printable
{
public:
  /// The base area's area ID
  std::string AreaID;
  /// The base area's parent ID
  std::string PrimaryAreaID;
  /// The base area's footprint
  Polygon Footprint;
  /// Properties within the base area
  std::vector<Property> Properties;
  /// Buildings within the base area
  std::vector<Building> Buildings;

  /// Pretty-print BaseArea.
  /// \return Pretty-print string
  std::string __str__() const override
  {
    return "Base area with area ID " + AreaID + ", primary area ID " +
           PrimaryAreaID + ", " + std::to_string(Properties.size()) +
           " properties and " + std::to_string(Buildings.size()) +
           " building(s)";
  }
};

} // namespace DTCCBUILDER

#endif // CORE_BASEAREA_H
