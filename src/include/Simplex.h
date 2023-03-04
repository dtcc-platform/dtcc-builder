// Simplex classes for 1D, 2D and 3D.
// Copyright (C) 2018 Anders Logg.
// Licensed under the MIT License

#ifndef DTCC_SIMPLEX_H
#define DTCC_SIMPLEX_H

#include <algorithm>
#include <tuple>
#include <vector>

namespace DTCC_BUILDER
{

class Simplex1D
{
public:
  // Vertex indices
  std::size_t v0{};
  std::size_t v1{};

  // Create default simplex
  Simplex1D() = default;

  // Create simplex with sorted vertices
  Simplex1D(std::size_t v0, std::size_t v1, bool sort = false)
  {
    /// Sort vertices if requested
    if (sort)
    {
      std::vector<size_t> v = {v0, v1};
      std::sort(v.begin(), v.end());
      this->v0 = v[0];
      this->v1 = v[1];
    }
    else
    {
      this->v0 = v0;
      this->v1 = v1;
    }
  }
};

class Simplex2D
{
public:
  // Vertex indices
  std::size_t v0{};
  std::size_t v1{};
  std::size_t v2{};

  // Create default simplex
  Simplex2D() = default;

  // Create simplex and optionally sort vertices
  Simplex2D(std::size_t v0, std::size_t v1, std::size_t v2,
            bool sort=false)
  {
    // Sort vertices if requested
    if (sort)
    {
      std::vector<size_t> v = {v0, v1, v2};
      std::sort(v.begin(), v.end());
      this->v0 = v[0];
      this->v1 = v[1];
      this->v2 = v[2];
    }
    else
    {
      this->v0 = v0;
      this->v1 = v1;
      this->v2 = v2;
    }
  }
};

class Simplex3D
{
public:
  // Vertex indices
  std::size_t v0{};
  std::size_t v1{};
  std::size_t v2{};
  std::size_t v3{};

  // Create default simplex
  Simplex3D() = default;

  // Create simplex and optionally sort vertices
  Simplex3D(std::size_t v0, std::size_t v1, std::size_t v2, std::size_t v3,
            bool sort=false)
  {
    // Sort vertices if requested
    if (sort)
    {
      std::vector<size_t> v = {v0, v1, v2, v3};
      std::sort(v.begin(), v.end());
      this->v0 = v[0];
      this->v1 = v[1];
      this->v2 = v[2];
      this->v3 = v[3];
    }
    else
    {
      this->v0 = v0;
      this->v1 = v1;
      this->v2 = v2;
      this->v3 = v3;
    }
  }
};

// Comparison function for Simplex1D (for use in e.g. maps)
struct CompareSimplex1D
{
  bool operator()(const Simplex1D &s, const Simplex1D &t) const
  {
    return std::tie(s.v0, s.v1) < std::tie(t.v0, t.v1);
  }
};

// Comparison function for Simplex2D (for use in e.g. maps)
struct CompareSimplex2D
{
  bool operator()(const Simplex2D &s, const Simplex2D &t) const
  {
    return std::tie(s.v0, s.v1, s.v2) < std::tie(t.v0, t.v1, t.v2);
  }
};

// Comparison function for Simplex3D (for use in e.g. maps)
struct CompareSimplex3D
{
  bool operator()(const Simplex3D &s, const Simplex3D &t) const
  {
    return std::tie(s.v0, s.v1, s.v2, s.v3) < std::tie(t.v0, t.v1, t.v2, t.v3);
  }
};

} // namespace DTCC_BUILDER

#endif
