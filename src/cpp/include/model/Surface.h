// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_SURFACE_H
#define DTCC_SURFACE_H

#include <vector>

#include "Logging.h"
#include "model/Vector.h"

namespace DTCC_BUILDER
{
class Surface : public Printable
{
public:
  std::vector<Vector3D> vertices{};
  std::vector<std::vector<Vector3D>> holes{};
  Surface() = default;
  virtual ~Surface() {} // make the destructor virtual

  /// Pretty-print
  std::string __str__() const override
  {
    return "Surface with " + str(vertices.size()) + " vertices";
  }
};

class MultiSurface
{
public:
  std::vector<Surface> surfaces{};

  MultiSurface() = default;
  virtual ~MultiSurface() {} // make the destructor virtual
};
} // namespace DTCC_BUILDER

#endif