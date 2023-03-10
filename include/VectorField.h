// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_VECTOR_FIELD_H
#define DTCC_VECTOR_FIELD_H

#include "Point.h"
#include "Vector.h"

namespace DTCC_BUILDER
{

  /// VectorField2D represents a vector field on a 2D domain.
  /// This is an interface class with several different implementations.
  class VectorField2D
  {
  public:

    /// Evaluate field at given point.
    ///
    /// @param p The point
    /// @return Value at point
    virtual Vector2D operator()(const Point2D& p) const = 0;

  };

  /// VectorField3D represents a vector field on a 3D domain.
  /// This is an interface class with several different implementations.
  class VectorField3D
  {
  public:

    /// Evaluate field at given point.
    ///
    /// @param p The point
    /// @return Value at point
    virtual Vector3D operator()(const Point3D& p) const = 0;

  };

}

#endif
