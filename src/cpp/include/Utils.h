// Copyright (C) 2020 Anders Logg, Anton J Olsson
// Licensed under the MIT License

#include <cassert>
#include <random>
#include <utility>

#ifndef DTCC_UTILS_H
#define DTCC_UTILS_H

namespace DTCC_BUILDER
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
  static std::string iso88591_to_utf8(std::basic_string<char> str)
  {
    std::string str_out;
    for (uint8_t ch : str)
    {
      if (ch < 0x80)
        str_out.push_back(ch);
      else
      {
        str_out.push_back(0xc0 | ch >> 6);
        str_out.push_back(0x80 | (ch & 0x3f));
      }
    }
    return str_out;
  }

  // Check if a string ends with another string
  static bool ends_with(const std::string &string, const std::string &ending)
  {
    if (ending.size() > string.size())
      return false;
    return std::equal(ending.rbegin(), ending.rend(), string.rbegin());
  }

  /// Return random number between 0 and 1
  static double random() { return std::rand() / double(RAND_MAX); }

  /// Returns the file name part of a path. If path end with a '/'
  /// return the directory name instead.
  /// \param str a string representing a path to a file.
  static std::string get_filename(std::string full_path,
                                  bool remove_extension = false)
  {
    if (full_path.back() == '/')
      full_path.pop_back();

    std::string file_name;
    size_t last_path_sep = full_path.find_last_of('/');

    if (last_path_sep != std::string::npos)
    {
      file_name = full_path.substr(last_path_sep + 1);
    }
    else // no directory
    {
      file_name = full_path;
    }
    if (remove_extension)
    {
      size_t file_ext_sep = file_name.find_last_of('.');
      if (file_ext_sep != std::string::npos)
      {
        file_name = file_name.substr(0, file_ext_sep);
      }
    }
    return file_name;
  }
};
} // namespace DTCC_BUILDER

#endif
