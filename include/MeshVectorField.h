// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_VECTOR_FIELD_H
#define DTCC_MESH_VECTOR_FIELD_H

#include "Mesh.h"
#include "VectorField.h"

namespace DTCC
{

  /// MeshVectorField2D represents a vector field on an irregular 2D mesh.
  /// The field can be efficiently evaluated at arbitrary points inside
  /// the mesh domain. The value is computed by linear interpolation.
  class MeshVectorField2D : public VectorField2D
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
    Point2D operator()(const Point2D& p) const
    {
      // FIXME: In progress
      return Point2D();
    }

  };

  /// MeshVectorField3D represents a vector field on an irregular 3D mesh.
  /// The field can be efficiently evaluated at arbitrary points inside
  /// the mesh domain. The value is computed by linear interpolation.
  class MeshVectorField3D : public VectorField3D
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
    Point3D operator()(const Point3D& p) const
    {
      // FIXME: In progress
      return Point3D();
    }

  };
}

#endif
