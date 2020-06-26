// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_COLOR_H
#define DTCC_COLOR_H

#include "Logging.h"

namespace DTCC
{

/// Color represents a color for BIM objects (RGB triplet).
class Color : public Printable
{
public:
  /// Red (R)
  size_t R{};

  /// Green (G)
  size_t G{};

  /// Blue (B)
  size_t B{};

  /// Empty constructor
  Color() {}

  /// Create color with given RGB values (0-255)
  Color(size_t R, size_t G, size_t B) : R(R), G(G), B(B) {}

  // Convert to string (pretty-print)
  std::string __str__() const
  {
    std::string s = "(" + str(R) + ", " + str(G) + ", " + str(B) + ")";
    return s;
  }
};

} // namespace DTCC

#endif
