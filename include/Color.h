// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_COLOR_H
#define DTCC_COLOR_H

#include <cstdint>

namespace DTCC
{
    class Color
    {
    public:
        uint8_t R{}; 
        uint8_t G{}; 
        uint8_t B{}; 
        uint8_t A{}; 

        Color(uint8_t r, uint8_t g,uint8_t b): R(r), G(g), B(b), A(255) {}
        Color(uint8_t r, uint8_t g,uint8_t b, uint8_t a): R(r), G(g), B(b), A(a) {}
    };

} // namespace DTCC
#endif
