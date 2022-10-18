// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_FIELD_H
#define DTCC_MESH_FIELD_H

#include "Mesh.h"
#include "Field.h"

namespace DTCC_BUILDER
{

  /// MeshField2D represents a scalar field on an irregular 2D mesh.
  /// The field can be efficiently evaluated at arbitrary points inside
  /// the mesh domain. The value is computed by linear interpolation.
  class MeshField2D : public Field2D
  {
  public:

    /// The mesh
    Mesh2D Mesh{};

    /// Array of values (vertex values)
    std::vector<double> Values{};

    /// Evaluate field at given point.
    ///
    /// @param p The point
    /// @return Value at point
    double operator()(const Vector2D& p) const
    {
      // FIXME: In progress
      return 0.0;
    }

  };

  /// MeshField3D represents a scalar field on an irregular 3D mesh.
  /// The field can be efficiently evaluated at arbitrary points inside
  /// the mesh domain. The value is computed by linear interpolation.
  class MeshField3D : public Field3D
  {
  public:

    /// The grid
    Mesh3D Mesh{};

    /// Array of values (vertex values)
    std::vector<double> Values{};

    /// Evaluate field at given point.
    ///
    /// @param p The point
    /// @return Value at point
    double operator()(const Vector3D& p) const
    {
      // FIXME: In progress
      return 0.0;
    }

  };
}

#endif
