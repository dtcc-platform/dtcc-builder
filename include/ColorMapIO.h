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
            
            // pixel 0 -> colormap value 1 and pixel height -> colormap value 0
            colorMap.InsertColor((num_colors-(i+1))/(num_colors-1.0), Color(r,g,b));     
        }
    }

    static void WritePNG(const ColorMap &colorMap, std::string fileName, size_t image_size = 256 ) {
        Color color;
        uint8_t r,g,b;
        png::image< png::rgb_pixel > image(image_size, image_size);
        for (size_t y = 0; y < image_size; y++)
        {
            double d = (image_size-(y+1))/(image_size-1.0);
            color = colorMap(d);
            r = (uint8_t) (color.R * 255);
            g = (uint8_t) (color.G * 255);
            b = (uint8_t) (color.B * 255);
            // std::cout << str(d) << ": "<< str(r) << ", " << str(r) <<std::endl;
            for (size_t x = 0; x < image_size; x++) {
                image.set_pixel(x,y,png::rgb_pixel(r,g,b));
            }
        }
        image.write(fileName); 
    }
};
}

#endif