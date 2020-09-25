// Copyright (C) 2020 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_COLOR_MAP_IO_FIELD_H
#define DTCC_COLOR_MAP_IO_FIELD_H

#include <png++/png.hpp>


#include "ColorMap.h"
#include "Logging.h"

namespace DTCC
{
class ColorMapIO
{
public:
    static void ReadPNG(ColorMap &colorMap, std::string fileName) {
        uint8_t r,g,b;        
        colorMap.Colors.clear();
        png::image< png::rgb_pixel > image(fileName);
        size_t num_colors = image.get_height();
        for (size_t i = 0; i < num_colors; i++)
        {
            auto p = image.get_pixel(0,i);
            r = (uint8_t) p.red;
            g = (uint8_t) p.green;
            b = (uint8_t) p.blue;
            colorMap.InsertColor(i/(num_colors-1.0), Color(r,g,b));     
        }
    }
};
}

#endif