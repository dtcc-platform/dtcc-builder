#ifndef DTCC_MESH_H
#define DTCC_MESH_H

#include <vector>

namespace DTCC_BUILDER
{

  class Point3D
  {
  public:
    double x{};
    double y{};
    double z{};
  };

  class Simplex3D
  {
  public:
    uint32_t v0{};
    uint32_t v1{};
    uint32_t v2{};
    uint32_t v3{};
  };

  class Mesh3D
  {
  public:
    std::vector<Point3D> Vertices{};
    std::vector<Simplex3D> Cells{};
  };

  class Mesh3DFlat
  {
  public:
    std::vector<double> Vertices{};
    std::vector<uint32_t> Cells{};
  };

}

#endif
