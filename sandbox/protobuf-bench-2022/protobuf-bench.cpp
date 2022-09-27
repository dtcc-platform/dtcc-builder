// Benchmarking and optimizing Protobuf set/get against DTCC native data types.
//
// To build without optimization:
//
//     cd build && cmake -DCMAKE_BUILD_TYPE=Debug ..
//
// To build with optimization:
//
//     cd build && cmake -DCMAKE_BUILD_TYPE=Release ..

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
void DTCCCreate(DTCC::Mesh3D &mesh, size_t N)
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

  //--- BENCH --------------------------------------
  mesh.Vertices.resize(numVertices);
  mesh.Cells.resize(numCells);
  //------------------------------------------------

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

	//--- BENCH --------------------------------------
	mesh.Vertices[k++] = DTCC::Point3D{x, y, z};
	//------------------------------------------------
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

	//--- BENCH --------------------------------------
	mesh.Cells[k++] = DTCC::Simplex3D{v0, v1, v3, v7};
	mesh.Cells[k++] = DTCC::Simplex3D{v0, v1, v7, v5};
	mesh.Cells[k++] = DTCC::Simplex3D{v0, v5, v7, v4};
	mesh.Cells[k++] = DTCC::Simplex3D{v0, v3, v2, v7};
	mesh.Cells[k++] = DTCC::Simplex3D{v0, v6, v4, v7};
	mesh.Cells[k++] = DTCC::Simplex3D{v0, v2, v6, v7};
	//------------------------------------------------
      }
    }
  }
}

// Compute volume of mesh
void DTCCAccess(const DTCC::Mesh3D &mesh)
{
  DTCC::Timer timer("DTCC access");

  double V{0.0};

  const auto &v{mesh.Vertices};

  for (const auto &c : mesh.Cells)
  {
    //--- BENCH --------------------------------------
    const auto &v0{v[c.v0]};
    const auto &v1{v[c.v1]};
    const auto &v2{v[c.v2]};
    const auto &v3{v[c.v3]};

    const double x0[]{v0.x, v0.y, v0.z};
    const double x1[]{v1.x, v1.y, v1.z};
    const double x2[]{v2.x, v2.y, v2.z};
    const double x3[]{v3.x, v3.y, v3.z};
    //------------------------------------------------

    V += TetrahedronVolume(x0, x1, x2, x3);
  }

  std::cout << "DTCC Volume = " << V << std::endl;
}

// Create N x N x N mesh of unit cube
void PROTCreate(PROT::Mesh3D &mesh, size_t N)
{
  DTCC::Timer timer("PROT create");

  const size_t nx = N;
  const size_t ny = N;
  const size_t nz = N;

  const double hx = 1.0 / static_cast<double>(nx);
  const double hy = 1.0 / static_cast<double>(ny);
  const double hz = 1.0 / static_cast<double>(nz);

  const size_t numVertices = (N+1)*(N+1)*(N+1);
  const size_t numCells = 6*N*N*N;

  //--- BENCH --------------------------------------
  //mesh.Vertices.resize(numVertices);
  //mesh.Cells.resize(numCells);
  //------------------------------------------------

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

	//--- BENCH --------------------------------------
	auto *vertex = mesh.add_vertices();
	vertex->set_x(x);
	vertex->set_y(y);
	vertex->set_z(z);
	//------------------------------------------------
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

	//--- BENCH --------------------------------------
	auto* c0 = mesh.add_cells();
	c0->set_v0(v0); c0->set_v1(v1); c0->set_v2(v3); c0->set_v3(v7);
	auto* c1 = mesh.add_cells();
	c1->set_v0(v0); c1->set_v1(v1); c1->set_v2(v7); c1->set_v3(v5);
	auto* c2 = mesh.add_cells();
	c2->set_v0(v0); c2->set_v1(v5); c2->set_v2(v7); c2->set_v3(v4);
	auto* c3 = mesh.add_cells();
	c3->set_v0(v0); c3->set_v1(v3); c3->set_v2(v2); c3->set_v3(v7);
	auto* c4 = mesh.add_cells();
	c4->set_v0(v0); c4->set_v1(v6); c4->set_v2(v4); c4->set_v3(v7);
	auto* c5 = mesh.add_cells();
	c5->set_v0(v0); c5->set_v1(v2); c5->set_v2(v6); c5->set_v3(v7);
	//------------------------------------------------
      }
    }
  }
}

// Compute volume of mesh
void PROTAccess(const PROT::Mesh3D &mesh)
{
  DTCC::Timer timer("PROT access");

  double V{0.0};

  for (size_t i = 0; i < mesh.cells_size(); i++)
  {
    //--- BENCH --------------------------------------
    const auto &c{mesh.cells(i)};
    const auto &v0{mesh.vertices(c.v0())};
    const auto &v1{mesh.vertices(c.v1())};
    const auto &v2{mesh.vertices(c.v2())};
    const auto &v3{mesh.vertices(c.v3())};

    const double x0[]{v0.x(), v0.y(), v0.z()};
    const double x1[]{v1.x(), v1.y(), v1.z()};
    const double x2[]{v2.x(), v2.y(), v2.z()};
    const double x3[]{v3.x(), v3.y(), v3.z()};
    //------------------------------------------------

    V += TetrahedronVolume(x0, x1, x2, x3);
  }

  std::cout << "PROT Volume = " << V << std::endl;
}

int main()
{
  // Size of mesh (unit cube N x N x N)
  const size_t N{256};

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
