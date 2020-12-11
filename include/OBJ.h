// Copyright (C) 2020 ReSpace AB
// Licensed under the MIT License

#ifndef DTCC_OBJ_H
#define DTCC_OBJ_H

#include <fstream>
#include <iostream>

#include "Parameters.h"
#include "Surface.h"

namespace DTCC
{

class OBJ
{
public:
  /// Write 3D surface to OBJ file.
  ///
  /// @param surface3D The 3D surface
  /// @parame fileName Filename (path)
  ///
  /// Note that the export flips the y and z coordinates:
  ///
  ///   -y --> z
  ///    z --> y
  ///
  /// so that the geometry may be correctly imported in Blender
  /// using standard Blender import settings (with -Z Forward, Y Up).
  ///
  /// Note also that Blender itself uses a standard mathematical coordinate
  /// system, but flips y and z coordinates by default during import.
  static void Write(const Surface3D &surface3D, std::string fileName)
  {
    Info("OBJ: Writing 3D surface file " + fileName + "...");

    // Open file
    std::ofstream f(fileName.c_str());
    if (!f)
    {
      Error("Unable to write to file: " + fileName);
    }

    // Write to string stream, much faster than writing directly to file
    std::stringstream s;

    // Set precision
    s << std::fixed << std::setprecision(Parameters::OutputPrecision);

    // Write header
    s << "# " << str(surface3D) << std::endl << std::endl;

    // Write vertices
    s << "# Vertices" << std::endl;
    for (size_t i = 0; i < surface3D.Vertices.size(); i++)
    {
      // Note flipped y and z coordinates!
      const Point3D &p = surface3D.Vertices[i];
      s << "v " << p.x << " " << p.z << " " << -p.y << std::endl;
    }
    s << std::endl;

    // Write normals
    s << "# Normals" << std::endl;
    for (size_t i = 0; i < surface3D.Vertices.size(); i++)
    {
      const Vector3D &n = surface3D.Normals[i];
      s << "vn " << n.x << " " << n.y << " " << n.z << std::endl;
    }
    s << std::endl;

    // Write triangles
    s << "# Triangles" << std::endl;
    for (size_t i = 0; i < surface3D.Faces.size(); i++)
    {
      const Simplex2D &f = surface3D.Faces[i];
      s << "f " << (f.v0 + 1) << " " << (f.v1 + 1) << " " << (f.v2 + 1)
        << std::endl;
    }

    // Write and close file
    f << s.str();
    f.close();
  }
};
}

#endif
