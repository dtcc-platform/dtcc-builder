// Copyright (C) 2020 Anders Logg, Anton J Olsson
// Licensed under the MIT License

#include <random>

#ifndef DTCC_UTILS_H
#define DTCC_UTILS_H

namespace DTCC
{

  /// Utils provides a collection of utility functions that are used
  /// throughout the library (string manipulation, command-line argument
  /// parsing, file management etc).
  class Utils
  {
  public:

    /// Crop integer x to interval [0, n - k). Requires care
    /// due to involving both signed and unsigned integers.
    ///
    /// @param x Signed integer
    /// @param n Unsigned integer (length of array)
    /// @param k Unsigned integer (maring at end of array)
    /// @return Unsigned integer within specified range
    static size_t crop(long int x, size_t n, size_t k = 0)
    {
      assert(n > 0);
      assert(k < n);
      return x < 0 ? 0 : (x + k >= n ? n - 1 - k : x);
    }

    /// Convert string encoded with ISO 8859-1 to UTF-8.
    /// \param str The ISO 8859-1 string to convert
    /// \return The string converted to UTF-8
    static std::string Iso88591ToUtf8(std::string &str)
    {
      std::string strOut;
      for (uint8_t ch : str)
      {
        if (ch < 0x80)
          strOut.push_back(ch);
        else
        {
          strOut.push_back(0xc0 | ch >> 6);
          strOut.push_back(0x80 | (ch & 0x3f));
        }
      }
      return strOut;
    }

    /// Return random number between 0 and 1
    static double Random() { return std::rand() / double(RAND_MAX); }
  };
}

#endif
