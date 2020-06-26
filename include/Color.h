// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_COLOR_H
#define DTCC_COLOR_H

#include "Logging.h"

#include <cstdint>

namespace DTCC
{

/// Colors are stored as doubles in the range 0-1
class Color
{
public:
  /// Red (R)
  double R{};

  /// Green (G)
  double G{};

  /// Blue (B)
  double B{};

  /// Alpha (A)
  double A{};

  /// Empty constructor
  Color() {}

  Color(double r, double g, double b) : R(r), G(g), B(b), A(1.0) {}
  Color(double r, double g, double b, double a) : R(r), G(g), B(b), A(a) {}

  // Convert to string (pretty-print)
  std::string __str__() const
  {
    std::string s = "(" + str(R) + ", " + str(G) + ", " + str(B) + ", " + str(A) +")";
    return s;
  }
};

} // namespace DTCC
#endif
