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

/* Does not seem to work for certain Protobuf version since parent class
   has been declared as final.

class Mesh3DProtFlat : public PROT::Mesh3DFlat
{
public:

  inline const double* Vertex(int32_t i) const
  {
    return &vertices()[3*i];
  }

};

*/

//--- DTCC NEST ----------------------------------------------------------------

void CreateDTCCNEST(DTCC::Mesh3D &mesh, uint32_t N)
{
  DTCC::Timer timer("Create DTCC NEST");

  const uint32_t nx = N;
  const uint32_t ny = N;
  const uint32_t nz = N;

  const double hx = 1.0 / static_cast<double>(nx);
  const double hy = 1.0 / static_cast<double>(ny);
  const double hz = 1.0 / static_cast<double>(nz);

  const uint32_t numVertices = (N+1)*(N+1)*(N+1);
  const uint32_t numCells = 6*N*N*N;

  //--- BENCH --------------------------------------
  mesh.Vertices.resize(numVertices);
  mesh.Cells.resize(numCells);
  //------------------------------------------------

  uint32_t k = 0;
  for (uint32_t iz = 0; iz < nz + 1; iz++)
  {
    for (uint32_t iy = 0; iy < ny + 1; iy++)
    {
      for (uint32_t ix = 0; ix < nx + 1; ix++)
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
  for (uint32_t iz = 0; iz < nz; iz++)
  {
    for (uint32_t iy = 0; iy < ny; iy++)
    {
      for (uint32_t ix = 0; ix < nx; ix++)
      {
	const uint32_t v0 = iz*(nx + 1)*(ny + 1) + iy*(nx + 1) + ix;
	const uint32_t v1 = v0 + 1;
	const uint32_t v2 = v0 + (nx + 1);
	const uint32_t v3 = v1 + (nx + 1);
	const uint32_t v4 = v0 + (nx + 1)*(ny + 1);
	const uint32_t v5 = v1 + (nx + 1)*(ny + 1);
	const uint32_t v6 = v2 + (nx + 1)*(ny + 1);
	const uint32_t v7 = v3 + (nx + 1)*(ny + 1);

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

void AccessDTCCNEST(const DTCC::Mesh3D &mesh)
{
  DTCC::Timer timer("Access DTCC NEST");

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

  std::cout << "DTCC NEST Volume = " << V << std::endl;
}

//--- DTCC FLAT ----------------------------------------------------------------

void CreateDTCCFLAT(DTCC::Mesh3DFlat &mesh, uint32_t N)
{
  DTCC::Timer timer("Create DTCC FLAT");

  const uint32_t nx = N;
  const uint32_t ny = N;
  const uint32_t nz = N;

  const double hx = 1.0 / static_cast<double>(nx);
  const double hy = 1.0 / static_cast<double>(ny);
  const double hz = 1.0 / static_cast<double>(nz);

  const uint32_t numVertices = (N+1)*(N+1)*(N+1);
  const uint32_t numCells = 6*N*N*N;

  //--- BENCH --------------------------------------
  mesh.Vertices.resize(3*numVertices);
  mesh.Cells.resize(4*numCells);
  //------------------------------------------------

  uint32_t k = 0;
  for (uint32_t iz = 0; iz < nz + 1; iz++)
  {
    for (uint32_t iy = 0; iy < ny + 1; iy++)
    {
      for (uint32_t ix = 0; ix < nx + 1; ix++)
      {
	const double x = ix*hx;
        const double y = iy*hy;
        const double z = iz*hz;

        //--- BENCH --------------------------------------
        mesh.Vertices[k++] = x;
        mesh.Vertices[k++] = y;
        mesh.Vertices[k++] = z;
        //------------------------------------------------
      }
    }
  }

  k = 0;
  for (uint32_t iz = 0; iz < nz; iz++)
  {
    for (uint32_t iy = 0; iy < ny; iy++)
    {
      for (uint32_t ix = 0; ix < nx; ix++)
      {
        const uint32_t v0 = iz*(nx + 1)*(ny + 1) + iy*(nx + 1) + ix;
        const uint32_t v1 = v0 + 1;
        const uint32_t v2 = v0 + (nx + 1);
        const uint32_t v3 = v1 + (nx + 1);
        const uint32_t v4 = v0 + (nx + 1)*(ny + 1);
        const uint32_t v5 = v1 + (nx + 1)*(ny + 1);
        const uint32_t v6 = v2 + (nx + 1)*(ny + 1);
        const uint32_t v7 = v3 + (nx + 1)*(ny + 1);

        //--- BENCH --------------------------------------
        mesh.Cells[k++] = v0; mesh.Cells[k++]= v1; mesh.Cells[k++] = v3; mesh.Cells[k++]= v7;
        mesh.Cells[k++] = v0; mesh.Cells[k++]= v1; mesh.Cells[k++] = v7; mesh.Cells[k++]= v5;
        mesh.Cells[k++] = v0; mesh.Cells[k++]= v5; mesh.Cells[k++] = v7; mesh.Cells[k++]= v4;
        mesh.Cells[k++] = v0; mesh.Cells[k++]= v3; mesh.Cells[k++] = v2; mesh.Cells[k++]= v7;
        mesh.Cells[k++] = v0; mesh.Cells[k++]= v6; mesh.Cells[k++] = v4; mesh.Cells[k++]= v7;
        mesh.Cells[k++] = v0; mesh.Cells[k++]= v2; mesh.Cells[k++] = v6; mesh.Cells[k++]= v7;
        //------------------------------------------------
      }
    }
  }
}

void AccessDTCCFLAT(const DTCC::Mesh3DFlat &mesh)
{
  DTCC::Timer timer("Access DTCC FLAT");

  double V{0.0};

  const auto &v{mesh.Vertices};
  const auto &c{mesh.Cells};

  const uint32_t n = mesh.Cells.size();
  for (uint32_t i = 0; i < n; i+= 4)
  {
    //--- BENCH --------------------------------------
    const uint32_t *vv = &c[i];
    const double *x0 = &v[3*vv[0]];
    const double *x1 = &v[3*vv[1]];
    const double *x2 = &v[3*vv[2]];
    const double *x3 = &v[3*vv[3]];

    //const auto v0{3*c[i]};
    //const auto v1{3*c[i + 1]};
    //const auto v2{3*c[i + 2]};
    //const auto v3{3*c[i + 3]};
    //const double x0[]{v[v0], v[v0 + 1], v[v0 + 2]};
    //const double x1[]{v[v1], v[v1 + 1], v[v1 + 2]};
    //const double x2[]{v[v2], v[v2 + 1], v[v2 + 2]};
    //const double x3[]{v[v3], v[v3 + 1], v[v3 + 2]};
    //------------------------------------------------

    V += TetrahedronVolume(x0, x1, x2, x3);
  }

  std::cout << "DTCC FLAT Volume = " << V << std::endl;
}

//--- PROT NEST ----------------------------------------------------------------

void CreatePROTNEST(PROT::Mesh3D &mesh, uint32_t N)
{
  DTCC::Timer timer("Create PROT NEST");

  const uint32_t nx = N;
  const uint32_t ny = N;
  const uint32_t nz = N;

  const double hx = 1.0 / static_cast<double>(nx);
  const double hy = 1.0 / static_cast<double>(ny);
  const double hz = 1.0 / static_cast<double>(nz);

  const uint32_t numVertices = (N+1)*(N+1)*(N+1);
  const uint32_t numCells = 6*N*N*N;

  //--- BENCH --------------------------------------
  //mesh.Vertices.resize(numVertices);
  //mesh.Cells.resize(numCells);
  //------------------------------------------------

  uint32_t k = 0;
  for (uint32_t iz = 0; iz < nz + 1; iz++)
  {
    for (uint32_t iy = 0; iy < ny + 1; iy++)
    {
      for (uint32_t ix = 0; ix < nx + 1; ix++)
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
  for (uint32_t iz = 0; iz < nz; iz++)
  {
    for (uint32_t iy = 0; iy < ny; iy++)
    {
      for (uint32_t ix = 0; ix < nx; ix++)
      {
	const uint32_t v0 = iz*(nx + 1)*(ny + 1) + iy*(nx + 1) + ix;
	const uint32_t v1 = v0 + 1;
	const uint32_t v2 = v0 + (nx + 1);
	const uint32_t v3 = v1 + (nx + 1);
	const uint32_t v4 = v0 + (nx + 1)*(ny + 1);
	const uint32_t v5 = v1 + (nx + 1)*(ny + 1);
	const uint32_t v6 = v2 + (nx + 1)*(ny + 1);
	const uint32_t v7 = v3 + (nx + 1)*(ny + 1);

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

void AccessPROTNEST(const PROT::Mesh3D &mesh)
{
  DTCC::Timer timer("Access PROT NEST");

  double V{0.0};

  const auto &v{mesh.vertices()};
  const auto &c{mesh.cells()};

  int32_t n = mesh.cells_size();
  for (int32_t i = 0; i < n; i++)
  {
    //--- BENCH --------------------------------------
    const auto &_c{c[i]};
    const auto &v0{v[_c.v0()]};
    const auto &v1{v[_c.v1()]};
    const auto &v2{v[_c.v2()]};
    const auto &v3{v[_c.v3()]};

    const double x0[]{v0.x(), v0.y(), v0.z()};
    const double x1[]{v1.x(), v1.y(), v1.z()};
    const double x2[]{v2.x(), v2.y(), v2.z()};
    const double x3[]{v3.x(), v3.y(), v3.z()};

    //const auto &_c{c[i]};
    //const auto &v0{v[_c.v0_]};
    //const auto &v1{v[_c.v1_]};
    //const auto &v2{v[_c.v2_]};
    //const auto &v3{v[_c.v3_]};

    //const double x0[]{v0.x_, v0.y_, v0.z_};
    //const double x1[]{v1.x_, v1.y_, v1.z_};
    //const double x2[]{v2.x_, v2.y_, v2.z_};
    //const double x3[]{v3.x_, v3.y_, v3.z_};
    //------------------------------------------------

    V += TetrahedronVolume(x0, x1, x2, x3);
  }

  std::cout << "PROT NEST Volume = " << V << std::endl;
}

//--- PROT FLAT  ---------------------------------------------------------------

void CreatePROTFLAT(PROT::Mesh3DFlat &mesh, uint32_t N)
{
  DTCC::Timer timer("Create PROT FLAT");

  const uint32_t nx = N;
  const uint32_t ny = N;
  const uint32_t nz = N;

  const double hx = 1.0 / static_cast<double>(nx);
  const double hy = 1.0 / static_cast<double>(ny);
  const double hz = 1.0 / static_cast<double>(nz);

  const uint32_t numVertices = (N+1)*(N+1)*(N+1);
  const uint32_t numCells = 6*N*N*N;

  //--- BENCH --------------------------------------
  //mesh.Vertices.resize(numVertices);
  //mesh.Cells.resize(numCells);
  //------------------------------------------------

  uint32_t k = 0;
  for (uint32_t iz = 0; iz < nz + 1; iz++)
  {
    for (uint32_t iy = 0; iy < ny + 1; iy++)
    {
      for (uint32_t ix = 0; ix < nx + 1; ix++)
      {
	const double x = ix*hx;
	const double y = iy*hy;
	const double z = iz*hz;

	//--- BENCH --------------------------------------
	mesh.add_vertices(x);
	mesh.add_vertices(y);
	mesh.add_vertices(z);
	//------------------------------------------------
      }
    }
  }

  k = 0;
  for (uint32_t iz = 0; iz < nz; iz++)
  {
    for (uint32_t iy = 0; iy < ny; iy++)
    {
      for (uint32_t ix = 0; ix < nx; ix++)
      {
	const uint32_t v0 = iz*(nx + 1)*(ny + 1) + iy*(nx + 1) + ix;
	const uint32_t v1 = v0 + 1;
	const uint32_t v2 = v0 + (nx + 1);
	const uint32_t v3 = v1 + (nx + 1);
	const uint32_t v4 = v0 + (nx + 1)*(ny + 1);
	const uint32_t v5 = v1 + (nx + 1)*(ny + 1);
	const uint32_t v6 = v2 + (nx + 1)*(ny + 1);
	const uint32_t v7 = v3 + (nx + 1)*(ny + 1);

	//--- BENCH --------------------------------------
	mesh.add_cells(v0); mesh.add_cells(v1); mesh.add_cells(v3); mesh.add_cells(v7);
	mesh.add_cells(v0); mesh.add_cells(v1); mesh.add_cells(v7); mesh.add_cells(v5);
	mesh.add_cells(v0); mesh.add_cells(v5); mesh.add_cells(v7); mesh.add_cells(v4);
	mesh.add_cells(v0); mesh.add_cells(v3); mesh.add_cells(v2); mesh.add_cells(v7);
	mesh.add_cells(v0); mesh.add_cells(v6); mesh.add_cells(v4); mesh.add_cells(v7);
        mesh.add_cells(v0); mesh.add_cells(v2); mesh.add_cells(v6); mesh.add_cells(v7);
	//------------------------------------------------
      }
    }
  }
}

void AccessPROTFLAT(PROT::Mesh3DFlat &mesh)
{
  DTCC::Timer timer("Access PROT FLAT");

  double V{0.0};

  const auto &v{mesh.vertices()};
  const auto &c{mesh.cells()};

  const int32_t n = mesh.cells_size();
  for (int32_t i = 0; i < n; i+= 4)
  {
    //--- BENCH --------------------------------------
    const uint32_t *vv = &c[i];
    //const double *x0 = mesh.Vertex(vv[0]);
    //const double *x1 = mesh.Vertex(vv[1]);
    //const double *x2 = mesh.Vertex(vv[2]);
    //const double *x3 = mesh.Vertex(vv[3]);
    const double *x0 = &v[3*vv[0]];
    const double *x1 = &v[3*vv[1]];
    const double *x2 = &v[3*vv[2]];
    const double *x3 = &v[3*vv[3]];

    //const auto v0{3*c[i]};
    //const auto v1{3*c[i + 1]};
    //const auto v2{3*c[i + 2]};
    //const auto v3{3*c[i + 3]};
    //const double x0[]{v[v0], v[v0 + 1], v[v0 + 2]};
    //const double x1[]{v[v1], v[v1 + 1], v[v1 + 2]};
    //const double x2[]{v[v2], v[v2 + 1], v[v2 + 2]};
    //const double x3[]{v[v3], v[v3 + 1], v[v3 + 2]};
    //------------------------------------------------

    V += TetrahedronVolume(x0, x1, x2, x3);
  }

  std::cout << "PROT FLAT Volume = " << V << std::endl;
}

int main()
{
  // Size of mesh (unit cube N x N x N)
  const uint32_t N{128};

  // DTCC NEST
  {
    DTCC::Mesh3D mesh{};
    CreateDTCCNEST(mesh, N);
    AccessDTCCNEST(mesh);
  }

  // DTCC FLAT
  {
    DTCC::Mesh3DFlat mesh{};
    CreateDTCCFLAT(mesh, N);
    AccessDTCCFLAT(mesh);
  }

  // PROT NEST
  {
    PROT::Mesh3D mesh{};
    CreatePROTNEST(mesh, N);
    AccessPROTNEST(mesh);
  }

  // PROT FLAT
  {
    PROT::Mesh3DFlat mesh{};
    CreatePROTFLAT(mesh, N);
    AccessPROTFLAT(mesh);
  }

  // Report timings
  DTCC::Timer::Report("protobuf-bench");

  return 0;
}
