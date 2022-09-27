#include <iostream>

#include "dtcc.pb.h"
#include "Mesh.h"
#include "Timer.h"

// Compute volume of tetrahedron
double TetrahedronVolume(const double x0[], const double x1[], const double x2[], const double x3[])
{
  return std::abs(x0[0]*(x1[1]*x2[2] + x3[1]*x1[2] + x2[1]*x3[2]
			 - x2[1]*x1[2] - x1[1]*x3[2] - x3[1]*x2[2])
		  - x1[0]*(x0[1]*x2[2] + x3[1]*x0[2] + x2[1]*x3[2]
			   - x2[1]*x0[2] - x0[1]*x3[2] - x3[1]*x2[2])
		  + x2[0]*(x0[1]*x1[2] + x3[1]*x0[2] + x1[1]*x3[2]
			   - x1[1]*x0[2] - x0[1]*x3[2] - x3[1]*x1[2]) -
		  x3[0]*(x0[1]*x1[2] + x1[1]*x2[2] + x2[1]*x0[2]
			 - x1[1]*x0[2] - x2[1]*x1[2] - x0[1]*x2[2])) / 6.0;
}

// Create N x N x N mesh of unit cube
void DTCCCreate(DTCC::Mesh3D& mesh, size_t N)
{
  DTCC::Timer timer("DTCC create");

  const size_t nx = N;
  const size_t ny = N;
  const size_t nz = N;

  const double hx = 1.0 / static_cast<double>(nx);
  const double hy = 1.0 / static_cast<double>(ny);
  const double hz = 1.0 / static_cast<double>(nz);

  const size_t numVertices = (N+1)*(N+1)*(N+1);
  const size_t numCells = 6*N*N*N;
  mesh.Vertices.resize(numVertices);
  mesh.Cells.resize(numCells);

  size_t k = 0;
  for (size_t iz = 0; iz < nz + 1; iz++)
  {
    for (size_t iy = 0; iy < ny + 1; iy++)
    {
      for (size_t ix = 0; ix < nx + 1; ix++)
      {
	const double x = ix*hx;
	const double y = iy*hy;
	const double z = iz*hz;
	mesh.Vertices[k++] = DTCC::Point3D{x, y, z};
      }
    }
  }

  k = 0;
  for (size_t iz = 0; iz < nz; iz++)
  {
    for (size_t iy = 0; iy < ny; iy++)
    {
      for (size_t ix = 0; ix < nx; ix++)
      {
	const size_t v0 = iz*(nx + 1)*(ny + 1) + iy*(nx + 1) + ix;
	const size_t v1 = v0 + 1;
	const size_t v2 = v0 + (nx + 1);
	const size_t v3 = v1 + (nx + 1);
	const size_t v4 = v0 + (nx + 1)*(ny + 1);
	const size_t v5 = v1 + (nx + 1)*(ny + 1);
	const size_t v6 = v2 + (nx + 1)*(ny + 1);
	const size_t v7 = v3 + (nx + 1)*(ny + 1);
	mesh.Cells[k++] = DTCC::Simplex3D{v0, v1, v3, v7};
	mesh.Cells[k++] = DTCC::Simplex3D{v0, v1, v7, v5};
	mesh.Cells[k++] = DTCC::Simplex3D{v0, v5, v7, v4};
	mesh.Cells[k++] = DTCC::Simplex3D{v0, v3, v2, v7};
	mesh.Cells[k++] = DTCC::Simplex3D{v0, v6, v4, v7};
	mesh.Cells[k++] = DTCC::Simplex3D{v0, v2, v6, v7};
      }
    }
  }
}

// Compute volume of mesh
void DTCCAccess(const DTCC::Mesh3D& mesh)
{
  DTCC::Timer timer("DTCC access");

  double x0[3]{}, x1[3]{}, x2[3]{}, x3[3]{};

  const auto & v{mesh.Vertices};
  double V{0.0};
  for (const auto &c : mesh.Cells)
  {
    const auto & v0{v[c.v0]};
    const auto & v1{v[c.v1]};
    const auto & v2{v[c.v2]};
    const auto & v3{v[c.v3]};

    const double x0[]{v0.x, v0.y, v0.z};
    const double x1[]{v1.x, v1.y, v1.z};
    const double x2[]{v2.x, v2.y, v2.z};
    const double x3[]{v3.x, v3.y, v3.z};
    V += TetrahedronVolume(x0, x1, x2, x3);
  }

  std::cout << "Volume = " << V << std::endl;
}

// Create N x N x N mesh of unit cube
void PROTCreate(PROT::Mesh3D& mesh, size_t N)
{
  DTCC::Timer timer("PROT create");
}

// Compute volume of mesh
void PROTAccess(const PROT::Mesh3D& mesh)
{
  DTCC::Timer timer("PROT access");
}

int main()
{
  // Size of mesh (unit cube N x N x N)
  const size_t N{128};

  // Run benchmarks for DTCC data structures
  {
    DTCC::Mesh3D mesh{};
    DTCCCreate(mesh, N);
    DTCCAccess(mesh);
  }

  // Run benchmarks for Protobuf data structures
  {
    PROT::Mesh3D mesh{};
    PROTCreate(mesh, N);
    PROTAccess(mesh);
  }

  // Report timings
  DTCC::Timer::Report("protobuf-bench");

  return 0;
}
