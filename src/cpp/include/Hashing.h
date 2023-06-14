// Copyright (C) 2020 ReSpace AB
// Licensed under the MIT License

#ifndef DTCC_HASHING_H
#define DTCC_HASHING_H

#include "model/Point.h"
#include "model/Simplices.h"
#include <iomanip>
#include <sstream>

namespace DTCC_BUILDER
{

class Hashing
{
public:
  /// Return hash of unsigned integer
  ///
  /// @param x Unsigned integer value
  /// @return Integer hash
  static size_t Hash(size_t x) { return std::hash<size_t>{}(x); }

  /// Return hash of double
  ///
  /// @param x Double value
  /// @return Integer hash
  static size_t Hash(double x) { return std::hash<double>{}(x); }

  /// Return hash of 2D simplex
  ///
  /// @param simplex Simplex2D
  /// @return Integer hash
  static size_t Hash(const Simplex2D &simplex)
  {
    const size_t h0 = Hashing::Hash(simplex.v0);
    const size_t h1 = Hashing::Hash(simplex.v1);
    const size_t h2 = Hashing::Hash(simplex.v2);
    return Hash(h0, h1, h2);
  }

  /// Return hash of 3D simplex
  ///
  /// @param simplex Simplex3D
  /// @return Integer hash
  static size_t Hash(const Simplex3D &simplex)
  {
    const size_t h0 = Hashing::Hash(simplex.v0);
    const size_t h1 = Hashing::Hash(simplex.v1);
    const size_t h2 = Hashing::Hash(simplex.v2);
    const size_t h3 = Hashing::Hash(simplex.v3);
    return Hash(h0, h1, h2, h3);
  }

  /// Return hash of 2D point
  ///
  /// @param point Point2D
  /// @return Integer hash
  static size_t Hash(const Point2D &point)
  {
    const size_t h0 = Hashing::Hash(point.x);
    const size_t h1 = Hashing::Hash(point.y);
    return Hash(h0, h1);
  }

  /// Return hash of 3D point
  ///
  /// @param point Point3D
  /// @return Integer hash
  static size_t Hash(const Point3D &point)
  {
    const size_t h0 = Hashing::Hash(point.x);
    const size_t h1 = Hashing::Hash(point.y);
    const size_t h2 = Hashing::Hash(point.z);
    return Hash(h0, h1, h2);
  }

  /// Combine 2 hashes
  ///
  /// @param h0 First hash
  /// @param h1 Second hash
  /// @return Combined hash
  static size_t Hash(size_t h0, size_t h1)
  {
    size_t seed = 0;
    HashCombine(seed, h0);
    HashCombine(seed, h1);
    return seed;
  }

  /// Combine 3 hashes
  ///
  /// @param h0 First hash
  /// @param h1 Second hash
  /// @param h2 Third hash
  /// @return Combined hash
  static size_t Hash(size_t h0, size_t h1, size_t h2)
  {
    size_t seed = 0;
    HashCombine(seed, h0);
    HashCombine(seed, h1);
    HashCombine(seed, h2);
    return seed;
  }

  /// Combine 4 hashes
  ///
  /// @param h0 First hash
  /// @param h1 Second hash
  /// @param h2 Third hash
  /// @param h3 Fourth hash
  /// @return Combined hash
  static size_t Hash(size_t h0, size_t h1, size_t h2, size_t h3)
  {
    size_t seed = 0;
    HashCombine(seed, h0);
    HashCombine(seed, h1);
    HashCombine(seed, h2);
    HashCombine(seed, h3);
    return seed;
  }

  /// Combine hashes (same as boost::hash_combine)
  static void HashCombine(size_t &seed, size_t hash)
  {
    seed ^= hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  /// Convert integer hash to hexadecimal string
  ///
  /// @param hash Integer hash (as returned by std::hash)
  /// @return String with hash in hexadecimal
  static std::string Hex(size_t hash)
  {
    std::stringstream s;
    s << std::setfill('0') << std::setw(sizeof(size_t) * 2) << std::hex << hash;
    return s.str();
  }
};

} // namespace DTCC_BUILDER

#endif
