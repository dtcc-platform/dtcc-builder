// Copyright (C) 2020 ReSpace AB
// Licensed under the MIT License

#ifndef DTCC_HASHING_H
#define DTCC_HASHING_H

#include "model/Simplices.h"
#include "model/Vector.h"
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
  static size_t hash(size_t x) { return std::hash<size_t>{}(x); }

  /// Return hash of double
  ///
  /// @param x Double value
  /// @return Integer hash
  static size_t hash(double x) { return std::hash<double>{}(x); }

  /// Return hash of 2D simplex
  ///
  /// @param simplex Simplex2D
  /// @return Integer hash
  static size_t hash(const Simplex2D &simplex)
  {
    const size_t h0 = Hashing::hash(simplex.v0);
    const size_t h1 = Hashing::hash(simplex.v1);
    const size_t h2 = Hashing::hash(simplex.v2);
    return hash(h0, h1, h2);
  }

  /// Return hash of 3D simplex
  ///
  /// @param simplex Simplex3D
  /// @return Integer hash
  static size_t hash(const Simplex3D &simplex)
  {
    const size_t h0 = Hashing::hash(simplex.v0);
    const size_t h1 = Hashing::hash(simplex.v1);
    const size_t h2 = Hashing::hash(simplex.v2);
    const size_t h3 = Hashing::hash(simplex.v3);
    return hash(h0, h1, h2, h3);
  }

  /// Return hash of 2D point
  ///
  /// @param point Vector2D
  /// @return Integer hash
  static size_t hash(const Vector2D &point)
  {
    const size_t h0 = Hashing::hash(point.x);
    const size_t h1 = Hashing::hash(point.y);
    return hash(h0, h1);
  }

  /// Return hash of 3D point
  ///
  /// @param point Vector3D
  /// @return Integer hash
  static size_t hash(const Vector3D &point)
  {
    const size_t h0 = Hashing::hash(point.x);
    const size_t h1 = Hashing::hash(point.y);
    const size_t h2 = Hashing::hash(point.z);
    return hash(h0, h1, h2);
  }

  /// Combine 2 hashes
  ///
  /// @param h0 First hash
  /// @param h1 Second hash
  /// @return Combined hash
  static size_t hash(size_t h0, size_t h1)
  {
    size_t seed = 0;
    hash_combine(seed, h0);
    hash_combine(seed, h1);
    return seed;
  }

  /// Combine 3 hashes
  ///
  /// @param h0 First hash
  /// @param h1 Second hash
  /// @param h2 Third hash
  /// @return Combined hash
  static size_t hash(size_t h0, size_t h1, size_t h2)
  {
    size_t seed = 0;
    hash_combine(seed, h0);
    hash_combine(seed, h1);
    hash_combine(seed, h2);
    return seed;
  }

  /// Combine 4 hashes
  ///
  /// @param h0 First hash
  /// @param h1 Second hash
  /// @param h2 Third hash
  /// @param h3 Fourth hash
  /// @return Combined hash
  static size_t hash(size_t h0, size_t h1, size_t h2, size_t h3)
  {
    size_t seed = 0;
    hash_combine(seed, h0);
    hash_combine(seed, h1);
    hash_combine(seed, h2);
    hash_combine(seed, h3);
    return seed;
  }

  /// Combine hashes (same as boost::hash_combine)
  static void hash_combine(size_t &seed, size_t hash)
  {
    seed ^= hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  /// Convert integer hash to hexadecimal string
  ///
  /// @param hash Integer hash (as returned by std::hash)
  /// @return String with hash in hexadecimal
  static std::string hex(size_t hash)
  {
    std::stringstream s;
    s << std::setfill('0') << std::setw(sizeof(size_t) * 2) << std::hex << hash;
    return s.str();
  }
};

} // namespace DTCC_BUILDER

#endif
