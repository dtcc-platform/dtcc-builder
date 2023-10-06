// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_BOUNDING_BOX_TREE_H
#define DTCC_BOUNDING_BOX_TREE_H

#include <limits>
#include <vector>

#include "BoundingBox.h"
#include "Geometry.h"
#include "Logging.h"

namespace DTCC_BUILDER
{

// Note: This is a new implementation of the corresponding algorithm
// in FEniCS/DOLFIN (by the same author) using simpler data structures
// and simpler design (separate 2D and 3D implementations).

/// BoundingBoxTree2D is a 2D axis-aligned bounding box tree (AABB tree).
class BoundingBoxTree2D : public Printable
{
public:
  // Tree node data structure. Each node has two children (unless it is a
  // leaf node) and a bounding box defining the region of the node.
  struct Node
  {
    // Index to first child node / -1 for leaf nodes
    int first{};

    // Index to second child node / object index for leaf nodes
    int second{};

    // Check if node is leaf
    bool is_leaf() const { return first == -1; }

    // Bounding box of node
    BoundingBox2D bbox;
  };

  // nodes of the bounding box tree. The tree is built bottom-up, meaning
  // that the leaf nodes will be added first and the top node last.
  std::vector<Node> nodes;

  // build bounding box tree for objects (defined by their bounding boxes)
  void build(const std::vector<BoundingBox2D> &bboxes)
  {
    // info("Building 2D bounding box tree for " +
    //     str(bboxes.size()) + " objects...");

    // clear tree if built before
    nodes.clear();

    // Skip if there is no data
    if (bboxes.empty())
      return;

    // Initialize indices of bounding boxes to be sorted
    std::vector<size_t> indices(bboxes.size());
    for (size_t i = 0; i < bboxes.size(); i++)
      indices[i] = i;

    // Recursively build bounding box tree
    build_recursive(bboxes, indices.begin(), indices.end());
  }

  /// find indices of bounding boxes containing given point.
  ///
  /// @param point The point
  /// @return Array of bounding box indices
  std::vector<size_t> find(const Vector2D &point) const
  {
    // Create empty list of bounding box indices
    std::vector<size_t> indices;

    // Recursively search bounding box tree for collisions
    find_recursive(indices, *this, point, nodes.size() - 1);

    return indices;
  }

  /// find indices of bounding boxes containing the 2d projection of a given
  /// point.
  ///
  /// @param point The point
  /// @return Array of bounding box indices
  std::vector<size_t> find(const Vector3D &point) const
  {
    return find(Vector2D(point.x, point.y));
  }

  /// find indices of bounding box collisions between this
  /// tree and given tree.
  ///
  /// @param tree The other tree
  /// @return Array of pairwise collisions
  std::vector<std::pair<size_t, size_t>>
  find(const BoundingBoxTree2D &tree) const
  {
    // Create empty list of bounding box indices
    std::vector<std::pair<size_t, size_t>> indices;

    // Recursively search bounding box tree for collisions
    find_recursive(indices, *this, tree, nodes.size() - 1,
                   tree.nodes.size() - 1);

    return indices;
  }

  /// Check if bounding box is empty
  bool empty() const { return nodes.empty(); }

  /// clear all data
  void clear() { nodes.clear(); }

  /// Pretty-print
  std::string __str__() const override
  {
    return "2D bounding box tree with " + str(nodes.size()) + " nodes";
  }

private:
  /// Comparison for partial sorting of bounding boxes along x-axis
  struct LessThanX
  {
    const std::vector<BoundingBox2D> &bboxes;
    explicit LessThanX(const std::vector<BoundingBox2D> &bboxes)
        : bboxes(bboxes)
    {
    }
    bool operator()(size_t i, size_t j)
    {
      return bboxes[i].P.x + bboxes[i].Q.x < bboxes[j].P.x + bboxes[j].Q.x;
    }
  };

  /// Comparison for partial sorting of bounding boxes along y-axis
  struct LessThanY
  {
    const std::vector<BoundingBox2D> &bboxes;
    explicit LessThanY(const std::vector<BoundingBox2D> &bboxes)
        : bboxes(bboxes)
    {
    }
    bool operator()(size_t i, size_t j)
    {
      return bboxes[i].P.y + bboxes[i].Q.y < bboxes[j].P.y + bboxes[j].Q.y;
    }
  };

  // build bounding box tree (recursive call). The input arguments are
  // the original full vector of bounding boxes and two iterators marking
  // the beginning and end of an array of indices (into the original vector
  // of bounding boxes) to be partitioned.
  int build_recursive(const std::vector<BoundingBox2D> &bboxes,
                      const std::vector<size_t>::iterator &begin,
                      const std::vector<size_t>::iterator &end)
  {
    // Create empty node
    Node node;

    // Check if we reached a leaf
    if (end - begin == 1)
    {
      node.bbox = bboxes[*begin];
      node.first = -1;
      node.second = *begin;
    }
    else
    {
      // Compute bounding box of bounding boxes
      node.bbox = compute_bounding_box(bboxes, begin, end);

      // Compute main axis of bounding box
      size_t main_axis = compute_main_axis(node.bbox);

      // Split boxes into two groups based on a partial sort along main axis
      auto middle = begin + (end - begin) / 2;
      if (main_axis == 0)
        std::nth_element(begin, middle, end, LessThanX(bboxes));
      else
        std::nth_element(begin, middle, end, LessThanY(bboxes));

      // Call recursively for both groups
      node.first = build_recursive(bboxes, begin, middle);
      node.second = build_recursive(bboxes, middle, end);
    }

    // Add node to tree
    nodes.push_back(node);

    // Return index of current node
    return nodes.size() - 1;
  }

  // find collisions between tree and point (recursive call)
  static void find_recursive(std::vector<size_t> &indices,
                             const BoundingBoxTree2D &tree,
                             const Vector2D &point,
                             size_t node_index)
  {
    // Get current node
    const Node &node = tree.nodes[node_index];

    // Check if node and point collide (if not, do nothing)
    if (Geometry::bounding_box_contains_2d(node.bbox, point))
    {
      // Check if leaf
      if (node.is_leaf())
      {
        // Add node index (we know that node collides with point)
        indices.push_back(node.second);
      }
      else
      {
        // Not a leaf node so check child nodes
        find_recursive(indices, tree, point, node.first);
        find_recursive(indices, tree, point, node.second);
      }
    }
  }

  // find collisions between tree A and tree B (recursive call)
  static void find_recursive(std::vector<std::pair<size_t, size_t>> &indices,
                             const BoundingBoxTree2D &tree_a,
                             const BoundingBoxTree2D &tree_b,
                             size_t node_index_a,
                             size_t node_index_b)
  {
    // Get current nodes
    const Node &node_a = tree_a.nodes[node_index_a];
    const Node &node_b = tree_b.nodes[node_index_b];

    // Check if nodes collide (if not, do nothing)
    if (Geometry::intersect_2d(node_a.bbox, node_b.bbox))
    {
      // Check if both nodes are leaves
      if (node_a.is_leaf() && node_b.is_leaf())
      {
        // Add node indices (we know that the nodes collide)
        indices.push_back(std::make_pair(node_a.second, node_b.second));
      }
      else if (node_a.is_leaf())
      {
        // If A is leaf, then descend B
        find_recursive(indices, tree_a, tree_b, node_index_a, node_b.first);
        find_recursive(indices, tree_a, tree_b, node_index_a, node_b.second);
      }
      else if (node_b.is_leaf())
      {
        // If B is leaf, then descend A
        find_recursive(indices, tree_a, tree_b, node_a.first, node_index_b);
        find_recursive(indices, tree_a, tree_b, node_a.second, node_index_b);
      }
      else
      {
        // If neither node is a leaf, traverse largest subtree
        if (node_index_a > node_index_b)
        {
          find_recursive(indices, tree_a, tree_b, node_a.first, node_index_b);
          find_recursive(indices, tree_a, tree_b, node_a.second, node_index_b);
        }
        else
        {
          find_recursive(indices, tree_a, tree_b, node_index_a, node_b.first);
          find_recursive(indices, tree_a, tree_b, node_index_a, node_b.second);
        }
      }
    }
  }

  // Compute bounding box of bounding boxes
  BoundingBox2D compute_bounding_box(const std::vector<BoundingBox2D> &bboxes,
                                     const std::vector<size_t>::iterator &begin,
                                     const std::vector<size_t>::iterator &end)
  {
    // Initialize bounding box
    constexpr double max = std::numeric_limits<double>::max();
    BoundingBox2D bounding_box;
    bounding_box.P = Vector2D(max, max);
    bounding_box.Q = Vector2D(-max, -max);

    // Iterate over bounding boxes to compute bounds
    for (auto it = begin; it != end; it++)
    {
      const BoundingBox2D &bbox = bboxes[*it];
      bounding_box.P.x = std::min(bounding_box.P.x, bbox.P.x);
      bounding_box.P.y = std::min(bounding_box.P.y, bbox.P.y);
      bounding_box.Q.x = std::max(bounding_box.Q.x, bbox.Q.x);
      bounding_box.Q.y = std::max(bounding_box.Q.y, bbox.Q.y);
    }

    return bounding_box;
  }

  // Compute main axis of bounding box
  size_t compute_main_axis(const BoundingBox2D &bbox)
  {
    const double dx = bbox.Q.x - bbox.P.x;
    const double dy = bbox.Q.y - bbox.P.y;
    return (dx > dy ? 0 : 1);
  }
};

/// BoundingBoxTree2D is a 2D axis-aligned bounding box tree (AABB tree).
class BoundingBoxTree3D
{
  // FIXME: Not yet implemented since not needed/used but can be easily
  // FIXME: implemented by extending the 2D version.
};

} // namespace DTCC_BUILDER

#endif
