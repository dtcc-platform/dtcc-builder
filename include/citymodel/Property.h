// Copyright (C) 2020 Anton J Olsson
// Licensed under the MIT License

#ifndef DTCC_PROPERTY_H
#define DTCC_PROPERTY_H

#include "Utils.h"
#include "citymodel/Building.h"
#include <utility>

namespace DTCC
{
/// Representation of a property.
class Property : public Printable
{
public:
  /// Property's UUID
  std::string UUID;

  /// Property's FNR
  size_t FNR;

  /// The property's total footprint
  Polygon Footprint;
  /// The buildings belonging to the property
  std::vector<Building> Buildings;

  /// Pretty-print Property.
  /// \return Pretty-print string.
  std::string __str__() const override
  {
    std::string prettyString = "Property with UUID " + UUID + ", FNR " +
                               str(FNR) + " and " + str(Buildings.size()) +
                               " building(s).\nFootprint: ";
    for (size_t i = 0; i < Footprint.Vertices.size(); ++i)
    {
      prettyString += Footprint.Vertices[i].__str__();
      if (i < Footprint.Vertices.size() - 1)
        prettyString += ", ";
    }
    return prettyString;
  }
};

} // namespace DTCC

#endif // DTCC_PROPERTY_H
