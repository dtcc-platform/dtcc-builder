/* Notes:
 *
 *  - NULL values will be removed. Are there to make partitions more intuitive.
 *  - Partition index "errors" are counting how many times index is getting out
 * of vertex column bounds.
 */

#ifndef DTCC_MESH_LAYERING_H
#define DTCC_MESH_LAYERING_H

#include "Logging.h"
#include "model/Mesh.h"

namespace DTCC_BUILDER
{

void color_sort(std::array<size_t, 3> &vertices, std::array<int, 3> &colors)
{
  // Combine colors and vertices into a vector of pairs
  std::vector<std::pair<size_t, int>> combinedList;
  for (size_t i = 0; i < colors.size(); ++i)
  {
    combinedList.emplace_back(vertices[i], colors[i]);
  }

  // Sort the combined list first based on the first element (vertices)
  std::sort(combinedList.begin(), combinedList.end(),
            [](const auto &lhs, const auto &rhs)
            { return lhs.first < rhs.first; });

  // Sort again based on the second element (colors)
  std::sort(combinedList.begin(), combinedList.end(),
            [](const auto &lhs, const auto &rhs)
            { return lhs.second < rhs.second; });

  // Unpack the sorted list back into colors and vertices
  for (size_t i = 0; i < combinedList.size(); ++i)
  {
    colors[i] = combinedList[i].second;
    vertices[i] = combinedList[i].first;
  }
}

VolumeMesh mesh_layering(const Mesh &mesh,
                         std::vector<double> &layer_heights,
                         std::vector<int> &face_colors,
                         const double domain_height)
{
  VolumeMesh volume_mesh;
  info("Mesh Layering Function.");
  const size_t num_vertices = mesh.vertices.size();
  const size_t num_faces = mesh.faces.size();

  const size_t num_layers = layer_heights.size();
  if (num_layers == 0)
  {
    error("Error: Empty Layer heights Vector.");
    return volume_mesh;
  }
  const double min_layer_height = layer_heights[0];
  const double max_layer_height = layer_heights[num_layers - 1];
  info("Number of layers: " + str(num_layers));
  info("Min Layer Height: " + str(min_layer_height));
  info("Max Layer Height: " + str(max_layer_height));

  const double adjusted_domain_height =
      std::ceil(domain_height / max_layer_height) * max_layer_height;
  info("Domain height adjusted to fit  chosen layer heights: " +
       str(adjusted_domain_height) + "m");

  // Each Vertex inherits the minimum color (layer height) from the faces it
  // belongs to.
  std::vector<int> vertex_colors(num_vertices, num_layers);
  for (size_t i = 0; i < num_faces; i++)
  {
    const std::vector<size_t> face = {mesh.faces[i].v0, mesh.faces[i].v1,
                                      mesh.faces[i].v2};
    for (const auto &j : face)
    {
      if (vertex_colors[j] > face_colors[i])
        vertex_colors[j] = face_colors[i];
    }
  }

  std::vector<std::vector<Vector3D>> col_vertices(num_vertices);
  std::vector<size_t> col_index_offset(num_vertices + 1, 0);

  for (size_t i = 0; i < num_vertices; i++)
  {
    const double layer_h = layer_heights[vertex_colors[i]];
    const size_t col_count =
        static_cast<size_t>((adjusted_domain_height / layer_h)) + 1;
    for (size_t j = 0; j < col_count; j++)
    {
      Vector3D v(mesh.vertices[i].x, mesh.vertices[i].y,
                 mesh.vertices[i].z + j * layer_h);
      col_vertices[i].push_back(v);
    }
    col_index_offset[i + 1] = col_index_offset[i] + col_vertices[i].size();
  }

  // for (size_t i = 0; i < 1000; i++)
  // {
  //   std::cout <<i <<") Col index offset "<< col_index_offset[i] <<" "
  //   <<col_vertices[i].size() <<std::endl;
  // }

  // Connecting vertices in each column to create cells for the volume mesh
  std::array<int, 4> partition_errors = {0};
  std::vector<Simplex3D> cells;
  for (size_t i = 0; i < num_faces; i++)
  {
    // std::cout << "Working on face "<< i+1 <<"/"<< num_faces<< std::endl;
    const Simplex2D face_simplex = mesh.faces[i];
    std::array<size_t, 3> face = {face_simplex.v0, face_simplex.v1,
                                  face_simplex.v2};
    std::array<int, 3> v_colors = {vertex_colors[face_simplex.v0],
                                   vertex_colors[face_simplex.v1],
                                   vertex_colors[face_simplex.v2]};
    color_sort(face, v_colors);

    const std::array<size_t, 3> column_offsets = {col_index_offset[face[0]],
                                                  col_index_offset[face[1]],
                                                  col_index_offset[face[2]]};
    // std::cout << "Face: "<<i << " coff0: "<< column_offsets[0] <<" coff1: "<<
    // column_offsets[1] << " coff2: "<< column_offsets[2] <<  std::endl;

    const std::array<size_t, 3> column_len = {
        col_index_offset[face[0] + 1] - col_index_offset[face[0]],
        col_index_offset[face[1] + 1] - col_index_offset[face[1]],
        col_index_offset[face[2] + 1] - col_index_offset[face[2]]};
    // std::cout << "Face: "<<i << " cl0: "<< column_len[0] <<" cl1: "<<
    // column_len[1] << " cl2: "<< column_len[2] <<  std::endl;
    const size_t num_prisms = column_len.back() - 1;

    // std::cout << "Face: "<<i << " Num of Prisms: "<< num_prisms << std::endl;
    std::vector<std::array<size_t, 3>> prism_iterator(num_prisms);

    for (size_t j = 0; j < num_prisms; j++)
    {
      prism_iterator[j] = {
          column_offsets[0] + j * (1 << (face_colors[i] - v_colors[0])),
          column_offsets[1] + j * (1 << (face_colors[i] - v_colors[1])),
          column_offsets[2] + j * (1 << (face_colors[i] - v_colors[2]))};
    }

    switch (3 * face_colors[i] - (v_colors[0] + v_colors[1] + v_colors[2]))
    {
    case 0:
    {
      for (const auto &ar : prism_iterator)
      {
        const size_t k = ar[0];
        const size_t l = ar[1];
        const size_t m = ar[2];

        if (k + 1 >= col_index_offset[face[0] + 1])
          partition_errors[0]++;
        if (l + 1 >= col_index_offset[face[1] + 1])
          partition_errors[0]++;
        if (m + 1 >= col_index_offset[face[2] + 1])
          partition_errors[0]++;

        std::array<size_t, 3> bot_triangle = {k, l, m};
        std::array<size_t, 3> mid_triangle = {NULL, NULL, NULL};
        std::array<size_t, 3> top_triangle = {k + 1, l + 1, m + 1};

        Simplex3D K0(bot_triangle[0], bot_triangle[1], bot_triangle[2],
                     top_triangle[2]);
        Simplex3D K1(bot_triangle[0], top_triangle[1], bot_triangle[1],
                     top_triangle[2]);
        Simplex3D K2(bot_triangle[0], top_triangle[0], top_triangle[1],
                     top_triangle[2]);

        cells.emplace_back(K0);
        cells.emplace_back(K1);
        cells.emplace_back(K2);
      }
    }
    break;
    case 1:
    {
      for (const auto &ar : prism_iterator)
      {
        const size_t k = ar[0];
        const size_t l = ar[1];
        const size_t m = ar[2];

        if (k + 2 >= col_index_offset[face[0] + 1])
          partition_errors[1]++;
        if (l + 1 >= col_index_offset[face[1] + 1])
          partition_errors[1]++;
        if (m + 1 >= col_index_offset[face[2] + 1])
          partition_errors[1]++;

        std::array<size_t, 3> bot_triangle = {k, l, m};
        std::array<size_t, 3> mid_triangle = {k + 1, NULL, NULL};
        std::array<size_t, 3> top_triangle = {k + 2, l + 1, m + 1};

        // Simplex3D K0(k, l, m, k + 1);
        // Simplex3D K1(l, m, k + 1, k + 1);
        // Simplex3D K2(l, m, k + 1, k + 2);
        // Simplex3D K3(k + 1, m + 1, k + 2, l + 1);

        Simplex3D K0(bot_triangle[0], bot_triangle[1], bot_triangle[2],
                     mid_triangle[0]);
        Simplex3D K1(bot_triangle[1], top_triangle[2], bot_triangle[2],
                     mid_triangle[0]);
        Simplex3D K2(bot_triangle[1], top_triangle[2], mid_triangle[0],
                     top_triangle[1]);
        Simplex3D K3(mid_triangle[0], top_triangle[0], top_triangle[1],
                     top_triangle[2]);

        cells.emplace_back(K0);
        cells.emplace_back(K1);
        cells.emplace_back(K2);
        cells.emplace_back(K3);
      }
    }
    break;
    case 2:
    {
      for (const auto &ar : prism_iterator)
      {
        const size_t k = ar[0];
        const size_t l = ar[1];
        const size_t m = ar[2];

        if (k + 2 >= col_index_offset[face[0] + 1])
        {
          partition_errors[2]++;
        }
        if (l + 2 >= col_index_offset[face[1] + 1])
        {
          partition_errors[2]++;
        }
        if (m + 1 >= col_index_offset[face[2] + 1])
        {
          partition_errors[2]++;
        }

        std::array<size_t, 3> bot_triangle = {k, l, m};
        std::array<size_t, 3> mid_triangle = {k + 1, l + 1, NULL};
        std::array<size_t, 3> top_triangle = {k + 2, l + 2, m + 1};

        // Simplex3D K0(k, l, m, l + 1);
        // Simplex3D K1(k, m, k + 1, l + 1);
        // Simplex3D K2(m, l, k + 1, l + 1);
        // Simplex3D K3(k + 2, l + 2, m + 1, k + 1);
        // Simplex3D K4(l + 2, m + 1, l + 1, k + 1);

        Simplex3D K0(bot_triangle[0], bot_triangle[1], bot_triangle[2],
                     mid_triangle[1]);
        Simplex3D K1(bot_triangle[0], bot_triangle[2], mid_triangle[0],
                     mid_triangle[1]);
        Simplex3D K2(bot_triangle[2], top_triangle[2], mid_triangle[0],
                     mid_triangle[1]);
        Simplex3D K3(top_triangle[0], top_triangle[2], top_triangle[1],
                     mid_triangle[0]);
        Simplex3D K4(top_triangle[1], top_triangle[2], mid_triangle[1],
                     mid_triangle[0]);

        cells.emplace_back(K0);
        cells.emplace_back(K1);
        cells.emplace_back(K2);
        cells.emplace_back(K3);
        cells.emplace_back(K4);
      }
    }
    break;

    case 3:
    {
      for (const auto &ar : prism_iterator)
      {
        const size_t k = ar[0];
        const size_t l = ar[1];
        const size_t m = ar[2];

        if (k + 2 >= col_index_offset[face[0] + 1])
        {
          // std::cout <<i << " Partition 3 k+2 = " << k+2 << " Next Column
          // index" <<  col_index_offset[face[0] + 1] << " | " <<
          // col_index_offset[face[0]] << std::endl;
          partition_errors[3]++;
          continue;
        }
        if (l + 2 >= col_index_offset[face[1] + 1])
        {
          // std::cout << "Partition 3 l+2 = " << l+2 << " Next Column index" <<
          // col_index_offset[face[1] + 1] << " | " << col_index_offset[face[1]]
          // << std::endl;
          partition_errors[3]++;
          continue;
        }
        if (m + 2 >= col_index_offset[face[2] + 1])
        {
          // std::cout << "Partition 3 m+2 = " << m+2 << " Next Column index" <<
          // col_index_offset[face[2] + 1] << " | " << col_index_offset[face[2]]
          // << std::endl;
          partition_errors[3]++;
          continue;
        }

        std::array<size_t, 3> bot_triangle = {k, l, m};
        std::array<size_t, 3> mid_triangle = {k + 1, l + 1, m + 1};
        std::array<size_t, 3> top_triangle = {k + 2, l + 2, m + 2};

        // Simplex3D K0(k, l, m, m + 1);
        // Simplex3D K1(k, m + 1, l, m + 1);
        // Simplex3D K2(k, m + 1, m, m + 1);
        // Simplex3D K3(m + 1, m + 1, m + 2, m + 2);
        // Simplex3D K4(m + 1, m + 2, m + 1, m + 2);
        // Simplex3D K5(m + 1, m + 2, m, m + 2);

        Simplex3D K0(bot_triangle[0], bot_triangle[1], bot_triangle[2],
                     mid_triangle[2]);
        Simplex3D K1(bot_triangle[0], mid_triangle[1], bot_triangle[1],
                     mid_triangle[2]);
        Simplex3D K2(bot_triangle[0], mid_triangle[0], mid_triangle[1],
                     mid_triangle[2]);
        Simplex3D K3(mid_triangle[0], mid_triangle[1], mid_triangle[2],
                     top_triangle[2]);
        Simplex3D K4(mid_triangle[0], top_triangle[1], mid_triangle[1],
                     top_triangle[2]);
        Simplex3D K5(mid_triangle[0], top_triangle[0], top_triangle[1],
                     top_triangle[2]);

        cells.emplace_back(K0);
        cells.emplace_back(K1);
        cells.emplace_back(K2);
        cells.emplace_back(K3);
        cells.emplace_back(K4);
        cells.emplace_back(K5);
      }
      break;
    }
    default:
      error("Face Coloring Error: Large layer height difference");
      break;
    }
  }
  std::cout << "Index skip error Cases:" << std::endl;
  std::cout << "Partition 0 Indexing errors: " << partition_errors[0]
            << std::endl;
  std::cout << "Partition 1 Indexing errors: " << partition_errors[1]
            << std::endl;
  std::cout << "Partition 2 Indexing errors: " << partition_errors[2]
            << std::endl;
  std::cout << "Partition 3 Indexing errors: " << partition_errors[3]
            << std::endl;

  std::vector<Vector3D> vertices;
  for (size_t j = 0; j < col_vertices.size(); j++)
  {
    for (size_t k = 0; k < col_vertices[j].size(); k++)
    {
      vertices.push_back(col_vertices[j][k]);
    }
  }

  volume_mesh.cells = cells;
  volume_mesh.vertices = vertices;
  return volume_mesh;
}

} // namespace DTCC_BUILDER

#endif // DTCC_MESH_LAYERING_H