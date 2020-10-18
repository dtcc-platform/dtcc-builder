// Copyright (C) 2020 Anton J Olsson
// Licensed under the MIT License

#ifndef DTCC_PROPERTY_H
#define DTCC_PROPERTY_H

#include "Building.h"
#include "Utils.h"
#include <utility>

namespace DTCC
{

class Property : public Printable
{
public:
  /// The property's total footprint.
  Polygon Footprint;
  /// The buildings belonging to the property.
  std::vector<Building> Buildings;

  /// Pretty-print.
  /// \return Pretty-print string.
  std::string __str__() const override
  {
    std::string prettyString = "Property with " +
                               std::to_string(Buildings.size()) +
                               " building(s).\n Footprint: ";
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
