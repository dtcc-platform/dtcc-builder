#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "CityBuilder.h"
#include "ElevationBuilder.h"
#include "MeshBuilder.h"
#include "MeshProcessor.h"
#include "MeshQualityMetrics.h"
#include "Smoother.h"
#include "VertexSmoother.h"
#include "model/Building.h"
#include "model/City.h"
#include "model/GridField.h"
#include "model/Mesh.h"
#include "model/PointCloud.h"
#include "model/Polygon.h"
#include "model/Simplices.h"
#include "model/Vector.h"
#include "Timer.h"

#include "include/FaceColoring.h"
#include "include/MeshLayering.h"

namespace py = pybind11;

namespace DTCC_BUILDER
{

VolumeMesh create_volume_mesh(py::array_t<double> vertices,
                 py::array_t<size_t> faces,
                 py::array_t<int> markers)
{
  Mesh mesh;
  auto verts_r = vertices.unchecked<2>();
  auto faces_r = faces.unchecked<2>();
  auto markers_r = markers.unchecked<1>();
  size_t num_vertices = verts_r.shape(0);
  size_t num_faces = faces_r.shape(0);
  size_t num_markers = markers_r.size();

  for (size_t i = 0; i < num_vertices; i++)
  {
    mesh.vertices.push_back(
        Vector3D(verts_r(i, 0), verts_r(i, 1), verts_r(i, 2)));
  }

  for (size_t i = 0; i < num_faces; i++)
  {
    mesh.faces.push_back(
        Simplex2D(faces_r(i, 0), faces_r(i, 1), faces_r(i, 2)));
  }

  for (size_t i = 0; i < num_markers; i++)
  {
    mesh.markers.push_back(markers_r(i));
  }
  info(mesh.__str__()) ;

  std::vector<double> layer_heights;
  std::vector<int> face_colors(num_faces) ; 

  Timer t1("Layer height computation");
  compute_layer_heights(mesh,layer_heights, face_colors);  
  t1.stop();
  t1.print();
  const double domain_height = 50; 

  Timer t2("Mesh layering");
  VolumeMesh vm = mesh_layering(mesh, layer_heights,face_colors, domain_height); 
  t2.stop();
  t2.print();
  return vm; 
}

} // namespace DTCC_BUILDER

PYBIND11_MODULE(_mesh_layering, m) {
    m.def("create_volume_mesh", &DTCC_BUILDER::create_volume_mesh, "Create C++ mesh");
}