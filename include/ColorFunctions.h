// Copyright (C) 2020 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_COLOR_FUNCTIONS_H
#define DTCC_COLOR_FUNCTIONS_H

#include <algorithm>
#include <vector>

#include "Color.h"
#include "Logging.h"
namespace DTCCBUILDER
{

class ColorFunctions
{
public:
  static Color Interpolate(Color c1, Color c2, double d)
  {
    double r,g,b;
    if (d <= 0)
      return c1;
    if (d >= 1)
      return c2;
    r = c1.R*(1-d)+c2.R*(d);
    g = c1.G*(1-d)+c2.G*(d);
    b = c1.B*(1-d)+c2.B*(d);
    return {r,g,b};
  }
};

} // namespace DTCCBUILDER
#endif