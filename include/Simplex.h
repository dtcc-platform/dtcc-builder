// Simplex classes for 2D and 3D.
// Copyright (C) 2018 Anders Logg.

#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <vector>

namespace VirtualCity
{

class Simplex1D
{
public:

  // vertex indices
  std::size_t v0;
  std::size_t v1;

  // Constructor
  Simplex1D(std::size_t v0,
            std::size_t v1)
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

  // vertex indices
  std::size_t v0;
  std::size_t v1;
  std::size_t v2;

  // Constructor
  Simplex2D(std::size_t v0,
            std::size_t v1,
            std::size_t v2)
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

  // vertex indices
  std::size_t v0;
  std::size_t v1;
  std::size_t v2;
  std::size_t v3;

  // Constructor
  Simplex3D(std::size_t v0,
            std::size_t v1,
            std::size_t v2,
            std::size_t v3)
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

}

#endif
