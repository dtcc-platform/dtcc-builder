// Copyright (C) 2020 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_COLOR_MAP_FIELD_H
#define DTCC_COLOR_MAP_FIELD_H

#include <vector>
#include <algorithm>

#include "Color.h"
#include "Logging.h"
namespace DTCC
{
class ColorMap
{
public:

    typedef std::pair<float,Color> colorMapEntry;
     
    std::vector<colorMapEntry> Colors{};

    ColorMap() {}

    void InsertColor(float startPoint, Color c)
    {
        Colors.push_back(std::make_pair(startPoint,c));
        sortColormap();
    }


    std::string __str__() const
    {
        std::string out = "Colormap: \n";
        for (auto c: Colors) {
            out += (str(c.first) + ": ");
            out += str(c.second);
            out += "\n";
        }
        return out;
    }

private:
    void sortColormap() {
        std::sort(std::begin(Colors), 
              std::end(Colors), 
              [](colorMapEntry a, colorMapEntry b ) {return a.first < b.first; });
    }
};
}
#endif