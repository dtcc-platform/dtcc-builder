// Copyright (C) 2020 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_COLOR_MAP_FIELD_H
#define DTCC_COLOR_MAP_FIELD_H

#include <vector>
#include <algorithm>
#include <stdexcept>


#include "Color.h"
#include "ColorFunctions.h"
#include "Logging.h"
namespace DTCC
{
enum ColorMapType {Linear,Discrete};

class ColorMap
{
public:

    typedef std::pair<float,Color> colorMapEntry;
    
    ColorMapType mapType;
    std::vector<colorMapEntry> Colors{};

    ColorMap(ColorMapType type = ColorMapType::Linear): mapType(type) {}

    void InsertColor(float startPoint, Color c)
    {
        if (startPoint<0 || startPoint>1)
        {
            throw std::invalid_argument("ColorMap data values must be between 0 and 1");
        }
        Colors.push_back(std::make_pair(startPoint,c));
        sortColormap();
    }

    size_t size() {
        return Colors.size();
    }

    Color operator()(double d) const
    {
        colorMapEntry lower;
        colorMapEntry higher;

        // if we are outside the range of the colormap 
        if (d<Colors.front().first)
            return Colors.front().second;
        if (d>Colors.back().first)
            return Colors.back().second;
        // inside the range of the color map 
        for (auto c: Colors) {
            if (d == c.first)
                return c.second;
            if (d > c.first)
                lower = c;
            if (d < c.first)
            {
                higher = c;
                break;
            }
        }
        if (mapType == Linear) {
            double interval_size = higher.first-lower.first;
            d = (d-lower.first)/interval_size;
            return ColorFunctions::Interpolate(lower.second,higher.second,d);
        } else
        {
            return lower.second;
        }
    }

    std::string __str__() const
    {
        std::string out = "Colormap: \n";
        for (auto c: Colors) 
        {
            out += (str(c.first) + ": ");   
            out += str(c.second);
            out += "\n";
        }
        return out;
    }

private:
    void sortColormap() 
    {
        std::sort(std::begin(Colors), 
              std::end(Colors), 
              [](colorMapEntry a, colorMapEntry b ) {return a.first < b.first; });
    }
};
}
#endif