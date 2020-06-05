// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

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

  };

}

#endif
