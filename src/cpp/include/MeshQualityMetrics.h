#ifndef DTCC_MESH_QUALITY_H
#define DTCC_MESH_QUALITY_H

#include <functional>
#include <iomanip>
#include <string>

#include "Logging.h"
#include "Table.h"
#include "model/Mesh.h"
#include "model/VolumeMesh.h"
//#include "QualityMetricsUtils.h"

namespace DTCC_BUILDER
{

// Constant Expressions used in different Metrics. (Could be moved to
// Constants.h)
static constexpr double pi_3 = 1.0471975512;
static constexpr double element_quality_tri = 6.92820323028; // = 4*sqrt(3)
static constexpr double element_quality_tetra = 124.70765802;
static constexpr double c_3_sqrt3_4 = 1.29903810568; // = 3*sqrt(3)/4

// Some Utility Functions used to compute or report the metrics.

inline std::string format_title(std::string str)
{
  std::string title_str = str + "\n";
  title_str += std::string(str.length(), '-') + "\n";
  return title_str;
}

inline void report(std::vector<double> &metric, std::string title_str)
{
  std::cout << format_title(title_str);

  std::string content_str = "";
  auto min = std::min_element(metric.begin(), metric.end());
  if (min != metric.end())
  {
    content_str += "MIN:\t" + std::to_string(*min) + "\t";
    // std::cout << "MIN:\t" << *min << std::endl;
  }
  auto max = std::max_element(metric.begin(), metric.end());
  if (max != metric.end())
  {
    content_str += "MAX:\t" + std::to_string(*max) + "\t";
    // std::cout << "MAX:\t" << *max << std::endl;
  }
  content_str += "\n" + std::string(content_str.length(), '-') + "\n";
  std::cout << content_str;
}

struct CellNeighbors
{
  std::vector<size_t> cells;
  std::vector<std::vector<size_t>> adjacent_faces;
};

template <typename T>
std::vector<T> vectorIntersection(const std::vector<T> &vector1,
                                  const std::vector<T> &vector2)
{
  std::vector<T> result;
  for (const T &element : vector1)
  {
    if (std::find(vector2.begin(), vector2.end(), element) != vector2.end())
    {
      result.push_back(element);
    }
  }
  return result;
}

std::vector<size_t> cell_intersection(Simplex3D s0, Simplex3D s1)
{
  std::vector<size_t> sv0 = {s0.v0, s0.v1, s0.v2, s0.v3};
  std::vector<size_t> sv1 = {s1.v0, s1.v1, s1.v2, s1.v3};
  return vectorIntersection(sv0, sv1);
}

std::vector<CellNeighbors> get_cell_adj_cells(const VolumeMesh &volume_mesh)
{
  std::cout << "Finding all neighbours of all Cells" << std::endl;
  const size_t num_cells = volume_mesh.cells.size();

  std::vector<CellNeighbors> neighbors(num_cells);
  for (size_t i = 0; i < num_cells; i++)
  {
    for (size_t j = i + 1; j < num_cells; j++)
    {
      if (neighbors[i].cells.size() == 4)
        break;
      std::vector<size_t> shared_vertices =
          cell_intersection(volume_mesh.cells[i], volume_mesh.cells[j]);

      if (shared_vertices.size() == 3)
      {
        neighbors[i].cells.push_back(j);
        neighbors[j].cells.push_back(i);

        neighbors[i].adjacent_faces.push_back(shared_vertices);
        neighbors[j].adjacent_faces.push_back(shared_vertices);
      }
    }
  }
  return neighbors;
}

// Class computing quality metrics for Triangular element (face)
class TriElementQuality
{
public:
  static double area(std::array<Vector3D, 3> &vertices)
  {
    Vector3D e0 = vertices[2] - vertices[1];
    Vector3D e1 = vertices[0] - vertices[2];
    Vector3D e2 = vertices[1] - vertices[0];
    const double l0 = e0.magnitude();
    const double l1 = e1.magnitude();
    const double l2 = e2.magnitude();

    double s = (l0 + l1 + l2) * 0.5;
    return sqrt(s * (s - l0) * (s - l1) * (s - l2));
  }

  static double element_quality(std::array<Vector3D, 3> &vertices)
  {

    double ar = area(vertices);

    Vector3D e0 = vertices[2] - vertices[1];
    Vector3D e1 = vertices[0] - vertices[2];
    Vector3D e2 = vertices[1] - vertices[0];

    const double l0 = e0.squared_magnitude();
    const double l1 = e1.squared_magnitude();
    const double l2 = e2.squared_magnitude();

    return element_quality_tri * ar / (l0 + l1 + l2);
  }

  static double edge_ratio(std::array<Vector3D, 3> &vertices)
  {
    Vector3D e0 = vertices[2] - vertices[1];
    Vector3D e1 = vertices[0] - vertices[2];
    Vector3D e2 = vertices[1] - vertices[0];

    double min_l = e0.magnitude();
    double max_l = e0.magnitude();

    if (max_l < e1.magnitude())
      max_l = e1.magnitude();
    if (max_l < e2.magnitude())
      max_l = e2.magnitude();

    if (min_l > e1.magnitude())
      min_l = e1.magnitude();
    if (min_l > e2.magnitude())
      min_l = e2.magnitude();

    if (min_l > 0)
      return max_l / min_l;
    return std::numeric_limits<double>::infinity();
  }

  static double aspect_ratio(std::array<Vector3D, 3> &vertices)
  {
    Vector3D e0 = vertices[2] - vertices[1];
    Vector3D e1 = vertices[0] - vertices[2];
    Vector3D e2 = vertices[1] - vertices[0];

    double l0 = e0.magnitude();
    double l1 = e1.magnitude();
    double l2 = e2.magnitude();

    double max_l = std::max({l0, l1, l2});
    double A = area(vertices);

    if (A > 0)
      return (max_l * (l0 + l1 + l2)) / (element_quality_tri * A);
    return std::numeric_limits<double>::infinity();
  }

  static double equiangular_skew(std::array<Vector3D, 3> &vertices)
  {
    std::array<double, 3> theta{};

    for (size_t i = 0; i < 3; i++)
    {
      Vector3D e0 = vertices[(i + 2) % 3] - vertices[i];
      Vector3D e1 = vertices[(i + 1) % 3] - vertices[i];
      theta[i] = e0.angle_between(e1);
    }

    auto theta_max = std::max_element(theta.begin(), theta.end());
    auto theta_min = std::min_element(theta.begin(), theta.end());

    double eq0 = (*theta_max - pi_3) / (2 * pi_3);
    double eq1 = (pi_3 - *theta_min) / pi_3;

    return std::max(eq0, eq1);
  }

  static double skewness(std::array<Vector3D, 3> &vertices)
  {
    Vector3D e0 = vertices[2] - vertices[1];
    Vector3D e1 = vertices[0] - vertices[2];
    Vector3D e2 = vertices[1] - vertices[0];

    double l0 = e0.magnitude();
    double l1 = e1.magnitude();
    double l2 = e2.magnitude();

    double A = area(vertices);

    double R = (l0 * l1 * l2) / (4 * A);
    double A_ideal = c_3_sqrt3_4 * R * R;

    return 1 - A / A_ideal;
  }
};

// Class containing quality metrics for surface meshes with triangular elements
class TriangleMeshQuality
{
private:
  struct Metric
  {
    std::function<std::vector<double>(const Mesh)> function;
    std::string title;
    Metric(std::function<std::vector<double>(const Mesh)> f, std::string t)
        : function(f), title(t){};
  };

  static std::vector<Metric> getAllMetricFunctions()
  {
    std::vector<Metric> metrics;
    metrics.push_back(Metric(area, "Face area (allowed to vary)"));
    metrics.push_back(
        Metric(element_quality, "Element quality (0-1, 1 is optimal)"));
    metrics.push_back(Metric(edge_ratio, "Edge ratio (1-inf, 1 is optimal)"));
    metrics.push_back(
        Metric(aspect_ratio, "Aspect ratio (1-inf, 1 is optimal)"));
    metrics.push_back(
        Metric(equiangular_skew, "Equiangular skewness (0-1, 0 is optimal)"));
    metrics.push_back(Metric(skewness, "Skewness (0-1, 0 is optimal)"));

    return metrics;
  };

public:
  static std::vector<double> area(const Mesh &mesh)
  {
    std::vector<double> area(mesh.faces.size(), 0);
    for (size_t f = 0; f < mesh.faces.size(); f++)
    {
      std::array<Vector3D, 3> vertices = {mesh.vertices[mesh.faces[f].v0],
                                          mesh.vertices[mesh.faces[f].v1],
                                          mesh.vertices[mesh.faces[f].v2]};
      area[f] = TriElementQuality::area(vertices);
    }
    return area;
  }

  static std::vector<double> element_quality(const Mesh &mesh)
  {
    std::vector<double> eq(mesh.faces.size(), 0);
    for (size_t f = 0; f < mesh.faces.size(); f++)
    {
      std::array<Vector3D, 3> vertices = {mesh.vertices[mesh.faces[f].v0],
                                          mesh.vertices[mesh.faces[f].v1],
                                          mesh.vertices[mesh.faces[f].v2]};
      eq[f] = TriElementQuality::element_quality(vertices);
    }
    return eq;
  }

  static std::vector<double> edge_ratio(const Mesh &mesh)
  {
    std::vector<double> er(mesh.faces.size(), 0);
    for (size_t f = 0; f < mesh.faces.size(); f++)
    {
      std::array<Vector3D, 3> vertices = {mesh.vertices[mesh.faces[f].v0],
                                          mesh.vertices[mesh.faces[f].v1],
                                          mesh.vertices[mesh.faces[f].v2]};
      er[f] = TriElementQuality::edge_ratio(vertices);
    }
    return er;
  }

  static std::vector<double> aspect_ratio(const Mesh &mesh)
  {
    std::vector<double> ar(mesh.faces.size(), 0);
    for (size_t f = 0; f < mesh.faces.size(); f++)
    {
      std::array<Vector3D, 3> vertices = {mesh.vertices[mesh.faces[f].v0],
                                          mesh.vertices[mesh.faces[f].v1],
                                          mesh.vertices[mesh.faces[f].v2]};
      ar[f] = TriElementQuality::aspect_ratio(vertices);
    }
    return ar;
  }

  static std::vector<double> equiangular_skew(const Mesh &mesh)
  {
    std::vector<double> sk(mesh.faces.size(), 0);
    for (size_t f = 0; f < mesh.faces.size(); f++)
    {
      std::array<Vector3D, 3> vertices = {mesh.vertices[mesh.faces[f].v0],
                                          mesh.vertices[mesh.faces[f].v1],
                                          mesh.vertices[mesh.faces[f].v2]};
      sk[f] = TriElementQuality::equiangular_skew(vertices);
    }
    return sk;
  }

  static std::vector<double> skewness(const Mesh &mesh)
  {
    std::vector<double> sk(mesh.faces.size(), 0);
    for (size_t f = 0; f < mesh.faces.size(); f++)
    {
      std::array<Vector3D, 3> vertices = {mesh.vertices[mesh.faces[f].v0],
                                          mesh.vertices[mesh.faces[f].v1],
                                          mesh.vertices[mesh.faces[f].v2]};
      sk[f] = TriElementQuality::skewness(vertices);
    }
    return sk;
  }

  static void quick_check_mesh(const Mesh &mesh)
  {
    {
      auto ar = aspect_ratio(mesh);
      auto max = std::max_element(ar.begin(), ar.end());
      if (max != ar.end())
        info("Max aspect ratio (1-inf, 1 is optimal):\t" +
             std::to_string(*max));
    }
    {
      auto sk = skewness(mesh);
      auto max = std::max_element(sk.begin(), sk.end());
      if (max != sk.end())
        info("Max skewness (0-1, 0 is optimal):\t" + std::to_string(*max));
    }
  }

  static void check_mesh(const Mesh &mesh, const int metrics_code)
  {
    if (metrics_code == 0)
    {
      quick_check_mesh(mesh);
      return;
    }

    int bit = 0;
    Table metric_table("Quality metrics for surface mesh");
    std::vector<std::string> header = {"Metric", "min", "max"};
    metric_table.rows.push_back(header);

    auto metrics_vector = getAllMetricFunctions();
    for (const auto &metric : metrics_vector)
    {
      if (metrics_code & 1 << bit)
      {

        std::vector<std::string> row;
        row.push_back(metric.title);

        auto q = metric.function(mesh);
        auto min = std::min_element(q.begin(), q.end());
        if (min != q.end())
        {
          row.push_back(std::to_string(*min));
        }
        auto max = std::max_element(q.begin(), q.end());
        if (max != q.end())
        {
          row.push_back(std::to_string(*max));
        }
        metric_table.rows.push_back(row);
      }
      bit++;
    }

    std::cout << metric_table.__str__() << std::endl;
  }
};

// Class computing quality metrics for Tetrahedron element (Cell)
class TetraElementQuality
{
public:
  static double volume(std::array<Vector3D, 4> &vertices)
  {
    Vector3D e01 = vertices[1] - vertices[0];
    Vector3D e02 = vertices[2] - vertices[0];
    Vector3D e03 = vertices[3] - vertices[0];

    double volume = std::abs(e01.dot(e02.cross(e03))) / 6;
    return volume;
  }

  static double element_quality(std::array<Vector3D, 4> &vertices)
  {
    Vector3D e01 = vertices[1] - vertices[0];
    Vector3D e02 = vertices[2] - vertices[0];
    Vector3D e03 = vertices[3] - vertices[0];
    Vector3D e12 = vertices[2] - vertices[1];
    Vector3D e13 = vertices[3] - vertices[1];
    Vector3D e23 = vertices[3] - vertices[2];

    double V = std::abs(e01.dot(e02.cross(e03))) / 6.0;
    double D =
        sqrt(std::pow(e01.squared_magnitude() + e02.squared_magnitude() +
                          e03.squared_magnitude() + e12.squared_magnitude() +
                          e13.squared_magnitude() + e23.squared_magnitude(),
                      3));

    return element_quality_tetra * (V / D);
  }

  static double skewness(std::array<Vector3D, 4> &vertices)
  {
    const double V = volume(vertices);
    Vector3D e01 = vertices[1] - vertices[0];
    Vector3D e02 = vertices[2] - vertices[0];
    Vector3D e03 = vertices[3] - vertices[0];
    Vector3D e12 = vertices[2] - vertices[1];
    Vector3D e13 = vertices[3] - vertices[1];
    Vector3D e23 = vertices[3] - vertices[2];

    // Circumradius
    const double p1 = e01.magnitude() * e23.magnitude();
    const double p2 = e02.magnitude() * e13.magnitude();
    const double p3 = e03.magnitude() * e12.magnitude();
    const double prod_1 = p1 + p2 + p3;
    const double prod_2 = p1 - p2 + p3;
    const double prod_3 = p1 + p2 - p3;
    const double prod_4 = -p1 + p2 + p3;

    double R = sqrt(prod_1 * prod_2 * prod_3 * prod_4) / (24 * V);

    double V_ideal = (8 * sqrt(3) / 27) * std::pow(R, 3);

    return 1 - (V / V_ideal);
  }

  static double edge_ratio(std::array<Vector3D, 4> &vertices)
  {
    const std::array<Vector3D, 6> edges = {
        vertices[1] - vertices[0], vertices[2] - vertices[0],
        vertices[3] - vertices[0], vertices[2] - vertices[1],
        vertices[3] - vertices[1], vertices[3] - vertices[2]};
    double max_l = edges[0].magnitude();
    double min_l = max_l;
    for (size_t i = 1; i < 6; i++)
    {
      double l = edges[i].magnitude();
      if (l > max_l)
        max_l = l;
      if (l < min_l)
        min_l = l;
    }
    if (min_l > 0)
      return max_l / min_l;
    return std::numeric_limits<double>::infinity();
  }
};

// Class containing quality metrics for volume meshes with Tetrahedral elements
class TetrahedronMeshQuality
{
private:
  struct Metric
  {
    std::function<std::vector<double>(const VolumeMesh)> function;
    std::string title;
    Metric(std::function<std::vector<double>(const VolumeMesh)> f,
           std::string t)
        : function(f), title(t){};
  };

  static std::vector<Metric> get_cell_metric_functions()
  {
    std::vector<Metric> metrics;
    metrics.push_back(Metric(volume, "Cell volume (allowed to vary)"));
    metrics.push_back(
        Metric(element_quality, "Element quality (0-1, 1 is optimal)"));
    metrics.push_back(Metric(edge_ratio, "Edge ratio (1-inf, 1 is optimal)"));
    metrics.push_back(Metric(skewness, "Skewness (0-1, 0 is optimal)"));

    return metrics;
  };

  static std::vector<Metric> get_face_metric_functions()
  {
    std::vector<Metric> metrics;
    metrics.push_back(Metric(area, "Face area (allowed to vary)"));
    metrics.push_back(Metric(size_change, "Size change (0-inf, 1 is optimal)"));
    metrics.push_back(
        Metric(orthogonal_quality, "Orthogonal quality (0-1, 1 is optimal)"));
    metrics.push_back(
        Metric(face_edge_ratio, "Face edge ratio (1-inf, 1 optimal)"));
    metrics.push_back(
        Metric(face_aspect_ratio, "Face aspect ratio (1-inf, 1 optimal)"));
    metrics.push_back(
        Metric(face_skewness, "Face skewness (0-1, 0 is optimal)"));

    return metrics;
  };

public:
  static std::vector<double> volume(const VolumeMesh &volume_mesh)
  {
    std::vector<double> vol(volume_mesh.cells.size(), 0);
    for (size_t c = 0; c < volume_mesh.cells.size(); c++)
    {
      std::array<Vector3D, 4> vertices = {
          volume_mesh.vertices[volume_mesh.cells[c].v0],
          volume_mesh.vertices[volume_mesh.cells[c].v1],
          volume_mesh.vertices[volume_mesh.cells[c].v2],
          volume_mesh.vertices[volume_mesh.cells[c].v3]};
      vol[c] = TetraElementQuality::volume(vertices);
    }
    return vol;
  }

  static std::vector<double> element_quality(const VolumeMesh &volume_mesh)
  {
    std::vector<double> eq(volume_mesh.cells.size(), 0);
    for (size_t c = 0; c < volume_mesh.cells.size(); c++)
    {
      std::array<Vector3D, 4> vertices = {
          volume_mesh.vertices[volume_mesh.cells[c].v0],
          volume_mesh.vertices[volume_mesh.cells[c].v1],
          volume_mesh.vertices[volume_mesh.cells[c].v2],
          volume_mesh.vertices[volume_mesh.cells[c].v3]};
      eq[c] = TetraElementQuality::element_quality(vertices);
    }
    return eq;
  }

  static std::vector<double> skewness(const VolumeMesh &volume_mesh)
  {
    std::vector<double> sk(volume_mesh.cells.size(), 0);
    for (size_t c = 0; c < volume_mesh.cells.size(); c++)
    {
      std::array<Vector3D, 4> vertices = {
          volume_mesh.vertices[volume_mesh.cells[c].v0],
          volume_mesh.vertices[volume_mesh.cells[c].v1],
          volume_mesh.vertices[volume_mesh.cells[c].v2],
          volume_mesh.vertices[volume_mesh.cells[c].v3]};
      sk[c] = TetraElementQuality::skewness(vertices);
    }
    return sk;
  }

  static std::vector<double> edge_ratio(const VolumeMesh &volume_mesh)
  {
    std::vector<double> er(volume_mesh.cells.size(), 0);
    for (size_t c = 0; c < volume_mesh.cells.size(); c++)
    {
      std::array<Vector3D, 4> vertices = {
          volume_mesh.vertices[volume_mesh.cells[c].v0],
          volume_mesh.vertices[volume_mesh.cells[c].v1],
          volume_mesh.vertices[volume_mesh.cells[c].v2],
          volume_mesh.vertices[volume_mesh.cells[c].v3]};
      er[c] = TetraElementQuality::edge_ratio(vertices);
    }
    return er;
  }

  static std::vector<double> size_change(const VolumeMesh &volume_mesh)
  {
    std::vector<double> size;
    auto neighbors = get_cell_adj_cells(volume_mesh);

    auto vol = volume(volume_mesh);

    for (size_t c = 0; c < volume_mesh.cells.size(); c++)
    {
      for (size_t adj_cell : neighbors[c].cells)
      {
        size.push_back(vol[c] / vol[adj_cell]);
      }
    }
    return size;
  }

  static std::vector<double> orthogonal_quality(const VolumeMesh &volume_mesh)
  {

    std::vector<double> oq(volume_mesh.cells.size());
    auto neighbors = get_cell_adj_cells(volume_mesh);
    for (size_t c = 0; c < volume_mesh.cells.size(); c++)
    {
      std::array<Vector3D, 4> vertices = {
          volume_mesh.vertices[volume_mesh.cells[c].v0],
          volume_mesh.vertices[volume_mesh.cells[c].v1],
          volume_mesh.vertices[volume_mesh.cells[c].v2],
          volume_mesh.vertices[volume_mesh.cells[c].v3]};

      Vector3D centroid =
          (vertices[0] + vertices[1] + vertices[2] + vertices[3]) / 4;
      std::vector<double> orthogonality;
      for (size_t face = 0; face < 4; face++)
      {
        std::array<Vector3D, 3> face_vertices = {vertices[(face + 1) % 4],
                                                 vertices[(face + 2) % 4],
                                                 vertices[(face + 3) % 4]};

        Vector3D e0 = vertices[1] - vertices[0];
        Vector3D e1 = vertices[2] - vertices[0];

        Vector3D face_centroid =
            (face_vertices[0] + face_vertices[1] + face_vertices[2]) / 3;
        Vector3D A = e0.cross(e1);
        A = A / A.magnitude();

        Vector3D f = face_centroid - centroid;
        f = f / f.magnitude();

        orthogonality.push_back(std::abs(A.dot(f)));
      }

      for (size_t n = 0; n < neighbors[c].cells.size(); n++)
      {
        Vector3D adj_centroid =
            (volume_mesh.vertices[volume_mesh.cells[neighbors[c].cells[n]].v0] +
             volume_mesh.vertices[volume_mesh.cells[neighbors[c].cells[n]].v1] +
             volume_mesh.vertices[volume_mesh.cells[neighbors[c].cells[n]].v2] +
             volume_mesh
                 .vertices[volume_mesh.cells[neighbors[c].cells[n]].v3]) /
            4;

        Vector3D e0 = volume_mesh.vertices[neighbors[c].adjacent_faces[n][1]] -
                      volume_mesh.vertices[neighbors[c].adjacent_faces[n][0]];
        Vector3D e1 = volume_mesh.vertices[neighbors[c].adjacent_faces[n][2]] -
                      volume_mesh.vertices[neighbors[c].adjacent_faces[n][0]];

        Vector3D A = e0.cross(e1);
        A = A / A.magnitude();

        Vector3D C = adj_centroid - centroid;
        C = C / C.magnitude();

        orthogonality.push_back(std::abs(A.dot(C)));
      }
      double skewness = TetraElementQuality::skewness(vertices);
      auto min_orthogonality =
          std::min_element(orthogonality.begin(), orthogonality.end());

      oq[c] = std::min(*min_orthogonality, skewness);
    }
    return oq;
  }

  static std::vector<double> area(const VolumeMesh &volume_mesh)
  {
    std::vector<double> a(volume_mesh.cells.size() * 4, 0);

    for (size_t c = 0; c < volume_mesh.cells.size(); c++)
    {
      const std::array<Vector3D, 4> vertices = {
          volume_mesh.vertices[volume_mesh.cells[c].v0],
          volume_mesh.vertices[volume_mesh.cells[c].v1],
          volume_mesh.vertices[volume_mesh.cells[c].v2],
          volume_mesh.vertices[volume_mesh.cells[c].v3]};

      for (size_t face = 0; face < 4; face++)
      {
        std::array<Vector3D, 3> v = {vertices[(face + 1) % 4],
                                     vertices[(face + 2) % 4],
                                     vertices[(face + 3) % 4]};

        a[4 * c + face] = TriElementQuality::area(v);
      }
    }
    return a;
  }

  static std::vector<double> face_edge_ratio(const VolumeMesh &volume_mesh)
  {
    std::vector<double> er(volume_mesh.cells.size() * 4, 0);

    for (size_t c = 0; c < volume_mesh.cells.size(); c++)
    {
      const std::array<Vector3D, 4> vertices = {
          volume_mesh.vertices[volume_mesh.cells[c].v0],
          volume_mesh.vertices[volume_mesh.cells[c].v1],
          volume_mesh.vertices[volume_mesh.cells[c].v2],
          volume_mesh.vertices[volume_mesh.cells[c].v3]};

      for (size_t face = 0; face < 4; face++)
      {
        std::array<Vector3D, 3> v = {vertices[(face + 1) % 4],
                                     vertices[(face + 2) % 4],
                                     vertices[(face + 3) % 4]};

        er[4 * c + face] = TriElementQuality::edge_ratio(v);
      }
    }
    return er;
  }

  static std::vector<double> face_aspect_ratio(const VolumeMesh &volume_mesh)
  {
    std::vector<double> ar(volume_mesh.cells.size() * 4, 0);

    for (size_t c = 0; c < volume_mesh.cells.size(); c++)
    {
      const std::array<Vector3D, 4> vertices = {
          volume_mesh.vertices[volume_mesh.cells[c].v0],
          volume_mesh.vertices[volume_mesh.cells[c].v1],
          volume_mesh.vertices[volume_mesh.cells[c].v2],
          volume_mesh.vertices[volume_mesh.cells[c].v3]};

      for (size_t face = 0; face < 4; face++)
      {
        std::array<Vector3D, 3> v = {vertices[(face + 1) % 4],
                                     vertices[(face + 2) % 4],
                                     vertices[(face + 3) % 4]};

        ar[4 * c + face] = TriElementQuality::aspect_ratio(v);
      }
    }
    return ar;
  }

  static std::vector<double> face_skewness(const VolumeMesh &volume_mesh)
  {
    std::vector<double> sk(volume_mesh.cells.size() * 4, 0);

    for (size_t c = 0; c < volume_mesh.cells.size(); c++)
    {
      const std::array<Vector3D, 4> vertices = {
          volume_mesh.vertices[volume_mesh.cells[c].v0],
          volume_mesh.vertices[volume_mesh.cells[c].v1],
          volume_mesh.vertices[volume_mesh.cells[c].v2],
          volume_mesh.vertices[volume_mesh.cells[c].v3]};

      for (size_t face = 0; face < 4; face++)
      {
        std::array<Vector3D, 3> v = {vertices[(face + 1) % 4],
                                     vertices[(face + 2) % 4],
                                     vertices[(face + 3) % 4]};

        sk[4 * c + face] = TriElementQuality::skewness(v);
      }
    }
    return sk;
  }

  static void quick_check_volume_mesh(const VolumeMesh &volume_mesh)
  {

    {
      auto ar = face_aspect_ratio(volume_mesh);
      auto max = std::max_element(ar.begin(), ar.end());
      if (max != ar.end())
        info("Face aspect ratio MAX (1-inf, 1 is optimal):\t" +
             std::to_string(*max));
    }
    {
      auto sk = skewness(volume_mesh);
      auto max = std::max_element(sk.begin(), sk.end());
      if (max != sk.end())
        info("Cell skewness MAX(0-1, 0 is optimal):\t" + std::to_string(*max));
    }
  }

  static void check_volume_mesh(const VolumeMesh &volume_mesh,
                                const int metrics_code)
  {
    if (metrics_code == 0)
    {
      quick_check_volume_mesh(volume_mesh);
      return;
    }

    int bit = 0;
    Table cell_metric_table("Cell-based quality metrics for volume mesh");
    std::vector<std::string> header = {"Metric", "min", "max"};
    cell_metric_table.rows.push_back(header);

    auto metrics_vector = get_cell_metric_functions();
    for (const auto &metric : metrics_vector)
    {

      if (metrics_code & 1 << bit)
      {
        std::vector<std::string> row;
        row.push_back(metric.title);
        auto q = metric.function(volume_mesh);
        auto min = std::min_element(q.begin(), q.end());
        if (min != q.end())
        {
          row.push_back(std::to_string(*min));
        }
        auto max = std::max_element(q.begin(), q.end());
        if (max != q.end())
        {
          row.push_back(std::to_string(*max));
        }
        cell_metric_table.rows.push_back(row);
      }
      bit++;
    }
    std::cout << cell_metric_table.__str__() << std::endl;

    Table face_metric_table("Face-based quality metrics for volume mesh");
    face_metric_table.rows.push_back(header);

    metrics_vector = get_face_metric_functions();
    for (const auto &metric : metrics_vector)
    {
      if (metrics_code & 1 << bit)
      {
        std::vector<std::string> row;
        row.push_back(metric.title);

        auto q = metric.function(volume_mesh);
        auto min = std::min_element(q.begin(), q.end());
        if (min != q.end())
        {
          row.push_back(std::to_string(*min));
        }
        auto max = std::max_element(q.begin(), q.end());
        if (max != q.end())
        {
          row.push_back(std::to_string(*max));
        }
        face_metric_table.rows.push_back(row);
      }
      bit++;
    }
    std::cout << face_metric_table.__str__() << std::endl;
  }
};

} // namespace DTCC_BUILDER
#endif