// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_FIELD_H
#define DTCC_FIELD_H

#include "Point.h"

namespace DTCC_BUILDER
{

  /// Field2D represents a scalar field on a 2D domain.
  /// This is an interface class with several different implementations.
  class Field2D
  {
  public:

    /// Evaluate field at given point.
    ///
    /// @param p The point
    /// @return Value at point
    virtual double operator()(const Point2D& p) const = 0;

    };

  /// Field3D represents a scalar field on a 3D domain.
  /// This is an interface class with several different implementations.
  class Field3D
  {
  public:

    /// Evaluate field at given point.
    ///
    /// @param p The point
    /// @return Value at point
    virtual double operator()(const Point3D& p) const = 0;

  };

}

#endif
