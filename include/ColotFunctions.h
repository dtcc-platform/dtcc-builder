// Copyright (C) 2020 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_COLOR_FUNCTIONS_H
#define DTCC_COLOR_FUNCTIONS_H

#include <algorithm>
#include <vector>

#include "Color.h"
#include "Logging.h"
namespace DTCC
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
    r = c1.R*d+c2.R*(1-d);
    g = c1.G*d+c2.G*(1-d);
    b = c1.B*d+c2.B*(1-d);
    return Color(r,g,b);
  }
};

} // namespace DTCC
#endif