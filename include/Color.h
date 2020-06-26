// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_COLOR_H
#define DTCC_COLOR_H

#include "Logging.h"

#include <cstdint>

namespace DTCC
{
class Color
{
public:
  /// Red (R)
  uint8_t R{};

  /// Green (G)
  uint8_t G{};

  /// Blue (B)
  uint8_t B{};

  /// Alpha (A)
  uint8_t A{};

  /// Empty constructor
  Color() {}

  Color(uint8_t r, uint8_t g, uint8_t b) : R(r), G(g), B(b), A(255) {}
  Color(uint8_t r, uint8_t g, uint8_t b, uint8_t a) : R(r), G(g), B(b), A(a) {}

  // Convert to string (pretty-print)
  std::string __str__() const
  {
    std::string s = "(" + str(R) + ", " + str(G) + ", " + str(B) + ", " + str(A) +")";
    return s;
  }
};

} // namespace DTCC
#endif
