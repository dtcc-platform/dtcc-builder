#ifndef DTCC_FACE_COLORING_H
#define DTCC_FACE_COLORING_H

#include "Logging.h"
#include "model/Mesh.h"

namespace DTCC_BUILDER
{

double ideal_layer_height(double area, const double scale = 1.0)
{
  // Compute ideal layer height for a regular tetrahedron
  const double c = std::pow(2.0, 1.5) * std::pow(3.0, -0.75);
  double h = scale * c * std::sqrt(area);
  return h;
}

double area2(std::array<Vector3D, 3> &vertices)
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

double area(std::array<Vector3D, 3> &vertices)
{
  const Vector3D e0 = vertices[2] - vertices[1];
  const Vector3D e1 = vertices[0] - vertices[2];
  const Vector3D e2 = vertices[0] - vertices[1];

  // Calculate the cross product of the two sides
  double crossProductX = e0.y * e2.z - e0.z * e2.y;
  double crossProductY = e0.z * e2.x - e0.x * e2.z;
  double crossProductZ = e0.x * e2.y - e0.y * e2.x;

  double area = 0.5 * std::sqrt(crossProductX * crossProductX +
                                crossProductY * crossProductY +
                                crossProductZ * crossProductZ);
  return area;
}

void assign_colors(std::vector<double> &layer_heights,
                   std::vector<double> &areas,
                   std::vector<int> &face_colors)
{
  const size_t num_faces = areas.size();

  // Assign layer heights to mesh (closest by quotient)
  for (size_t i = 0; i < num_faces; i++)
  {
    double h = ideal_layer_height(areas[i]);
    std::vector<double> d(layer_heights.size());
    d.reserve(layer_heights.size());

    for (size_t ih = 0; ih < layer_heights.size(); ih++)
    {
      d[ih] = std::abs(std::log(h / layer_heights[ih]));
    }
    auto min_it = std::min_element(d.begin(), d.end());
    face_colors[i] = std::distance(d.begin(), min_it);
  }
}

void reassign_colors(std::vector<int> &face_colors,
                     std::vector<std::unordered_set<int>> &ff)
{
  // Reassign colors to avoid big jumps
  for (size_t i = 0; i < ff.size(); i++)
  {
    for (const auto &j : ff[i])
    {
      int diff = face_colors[i] - face_colors[j];
      if (diff > 1)
      {
        face_colors[i] -= diff - 1;
      }
    }
  }
}

size_t check_neighbors(std::vector<int> &face_colors,
                       std::vector<std::unordered_set<int>> &ff)
{
  int max_diff = 0;
  size_t num_big_diffs = 0;
  std::vector<int> big_diff_colors;
  for (size_t i = 0; i < ff.size(); i++)
  {
    const auto face = ff[i];
    int diff_num = 0;
    for (const auto &f : face)
    {
      int diff = face_colors[i] - face_colors[f];
      if (diff > 1)
      {
        max_diff = std::max(max_diff, diff);
        diff_num++;
      }
    }
    if (diff_num > 0)
    {
      num_big_diffs++;
      big_diff_colors.push_back(1);
    }
    else
    {
      big_diff_colors.push_back(0);
    }
  }
  double big_diff_percentage =
      100.0 * num_big_diffs / static_cast<double>(ff.size());
  const std::string s = "Num Big diffs " + str(num_big_diffs) + " / " +
                        str(ff.size()) + " (" + str(big_diff_percentage, 2L) +
                        "%)";

  info(s);
  return num_big_diffs;
}

void compute_layer_heights(const Mesh &mesh,
                           std::vector<double> &layer_heights,
                           std::vector<int> &face_colors)
{
  std::vector<double> areas(mesh.faces.size());
  for (std::size_t i = 0; i < mesh.faces.size(); i++)
  {
    std::array<Vector3D, 3> vertices = {mesh.vertices[mesh.faces[i].v0],
                                        mesh.vertices[mesh.faces[i].v1],
                                        mesh.vertices[mesh.faces[i].v2]};
    areas[i] = area(vertices);
  }

  // Compute ideal layer heights for smallest and largest mesh sizes
  double _min_area = 0.0;
  auto min = std::min_element(areas.begin(), areas.end());
  if (min != areas.end())
  {
    _min_area = *min;
  }

  double _max_area = 0.0;
  auto max = std::max_element(areas.begin(), areas.end());
  if (max != areas.end())
  {
    _max_area = *max;
  }

  info("Min Face Area: " + str(_min_area));
  info("Max Face Area: " + str(_max_area));

  info("Min Face Ideal Layer Height: " + str(ideal_layer_height(_min_area)));
  info("Max Face Ideal Layer Height: " + str(ideal_layer_height(_max_area)));
  // info("Test Layer Height: " + str(ideal_layer_height(1.0)));

  double _min_height = ideal_layer_height(_min_area);
  double _max_height = ideal_layer_height(_max_area);

  // Compute dyadic mesh sizes to match min/max as close as possible
  double rho = _max_height / _min_height;
  double mid = std::sqrt(_min_height * _max_height);
  int num_layers = static_cast<int>(std::log2(rho) + 0.5);
  double min_height = mid / std::pow(2, num_layers / 2);

  info("rho: " + str(rho) + " | mid: " + str(mid));
  info("Num of Layers: " + str(num_layers));

  // Create layer_heights array
  for (int i = 0; i < num_layers; i++)
  {
    layer_heights.push_back(min_height * std::pow(2.0, i));
  }

  for (size_t i = 0; i < layer_heights.size(); i++)
  {
    std::cout << "Layer " << i << ": " << layer_heights[i] << "m" << std::endl;
  }

  info("Assign layer heights to mesh");
  // Assign layer heights to mesh (closest by quotient)
  assign_colors(layer_heights, areas, face_colors);

  info("Build mapping from vertices to faces");
  // Build mapping from vertices to faces
  std::vector<std::unordered_set<int>> vf(mesh.vertices.size());
  for (size_t i = 0; i < mesh.faces.size(); i++)
  {
    vf[mesh.faces[i].v0].insert(i);
    vf[mesh.faces[i].v1].insert(i);
    vf[mesh.faces[i].v2].insert(i);
  }

  info("Build mapping from faces to faces");
  // Build mapping from faces to faces
  std::vector<std::unordered_set<int>> ff(mesh.faces.size());
  for (size_t i = 0; i < mesh.faces.size(); i++)
  {
    for (size_t j = 0; j < 3; j++)
    {
      for (const int &v : vf[mesh.faces[i][j]])
      {
        ff[i].insert(v);
      }
    }
  }

  info("Iteratively reassign colors to avoid big jumps");
  // Iteratively reassign colors to avoid big jumps
  // Note: It should not do 3 iterations but reassign colors until there are
  // 0 big jumps.
  for (size_t i = 0; i < 3; i++)
  {
    if (check_neighbors(face_colors, ff))
    {
      info("Reassigning colors, iteration " + str(i));
      reassign_colors(face_colors, ff);
    }
  }
}

} // namespace DTCC_BUILDER

#endif // DTCC_FACE_COLORING_H