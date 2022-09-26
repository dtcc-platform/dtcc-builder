#ifndef DTCC_MESH_H
#define DTCC_MESH_H

#include <vector>

namespace DTCC
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
    std::size_t v0{};
    std::size_t v1{};
    std::size_t v2{};
    std::size_t v3{};
  };

  class Mesh3D
  {
  public:

    std::vector<Point3D> Vertices{};
    std::vector<Simplex3D> Cells{};
    std::vector<int> Markers{};

  };

}

#endif
