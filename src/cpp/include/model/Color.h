// Copyright (C) 2020 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_COLOR_H
#define DTCC_COLOR_H

#include "Logging.h"

#include <cstdint>

namespace DTCC_BUILDER
{

/// colors are stored as doubles in the range 0-1
class Color : public Printable
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

  /// empty constructor
  Color() = default;
  virtual ~Color() {} // make the destructor virtual

  Color(double r, double g, double b) : R(r), G(g), B(b), A(1.0) {}
  Color(double r, double g, double b, double a) : R(r), G(g), B(b), A(a) {}

  Color(uint8_t r, uint8_t g, uint8_t b)
      : R(r / 255.0), G(g / 255.0), B(b / 255.0), A(1.0)
  {
  }
  Color(uint8_t r, uint8_t g, uint8_t b, uint8_t a)
      : R(r / 255.0), G(g / 255.0), B(b / 255.0), A(a / 255.0)
  {
  }

  // Convert to string (pretty-print)
  std::string __str__() const override
  {
    std::string s =
        "(" + str(R) + ", " + str(G) + ", " + str(B) + ", " + str(A) + ")";
    return s;
  }
};

} // namespace DTCC_BUILDER
#endif
