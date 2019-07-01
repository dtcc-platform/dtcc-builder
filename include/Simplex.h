// Simplex classes for 1D, 2D and 3D.
// Copyright (C) 2018 Anders Logg.

#ifndef VC_SIMPLEX_H
#define VC_SIMPLEX_H

#include <vector>

namespace VirtualCity
{

class Simplex1D
{
public:
  // Vertex indices
  std::size_t v0 {};
  std::size_t v1 {};

  // Create default simplex
  Simplex1D() {};

  // Create simplex with sorted vertices
  Simplex1D(std::size_t v0, std::size_t v1)
  {
    /// Sort vertices
    std::vector<size_t> v = {v0, v1};
    std::sort(v.begin(), v.end());
    this->v0 = v[0];
    this->v1 = v[1];
  }
};

class Simplex2D
{
public:
  // Vertex indices
  std::size_t v0 {};
  std::size_t v1 {};
  std::size_t v2 {};

  // Create default simplex
  Simplex2D() {};

  // Create simplex with sorted vertices
  Simplex2D(std::size_t v0, std::size_t v1, std::size_t v2)
  {
    // Sort vertices
    std::vector<size_t> v = {v0, v1, v2};
    std::sort(v.begin(), v.end());
    this->v0 = v[0];
    this->v1 = v[1];
    this->v2 = v[2];
  }
};

class Simplex3D
{
public:
  // Vertex indices
  std::size_t v0 {};
  std::size_t v1 {};
  std::size_t v2 {};
  std::size_t v3 {};

  // Create default simplex
  Simplex3D() {};

  // Create simplex with sorted vertices
  Simplex3D(std::size_t v0, std::size_t v1, std::size_t v2, std::size_t v3)
  {
    // Sort vertices
    std::vector<size_t> v = {v0, v1, v2, v3};
    std::sort(v.begin(), v.end());
    this->v0 = v[0];
    this->v1 = v[1];
    this->v2 = v[2];
    this->v3 = v[3];
  }
};

} // namespace VirtualCity

#endif
