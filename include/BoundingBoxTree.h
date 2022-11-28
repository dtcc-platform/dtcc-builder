// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_BOUNDING_BOX_TREE_H
#define DTCC_BOUNDING_BOX_TREE_H

#include <vector>
#include <limits>

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
    bool IsLeaf() const { return first == -1; }

    // Bounding box of node
    BoundingBox2D bbox;
  };

  // Nodes of the bounding box tree. The tree is built bottom-up, meaning
  // that the leaf nodes will be added first and the top node last.
  std::vector<Node> Nodes;

  // Build bounding box tree for objects (defined by their bounding boxes)
  void Build(const std::vector<BoundingBox2D> &bboxes)
  {
    // info("BoundingBoxTree: Building 2D bounding box tree for " +
    //     str(bboxes.size()) + " objects...");

    // Clear tree if built before
    Nodes.clear();

    // Skip if there is no data
    if (bboxes.empty())
      return;

    // Initialize indices of bounding boxes to be sorted
    std::vector<size_t> indices(bboxes.size());
    for (size_t i = 0; i < bboxes.size(); i++)
      indices[i] = i;

    // Recursively build bounding box tree
    BuildRecursive(bboxes, indices.begin(), indices.end());
  }

  /// Find indices of bounding boxes containing given point.
  ///
  /// @param point The point
  /// @return Array of bounding box indices
  std::vector<size_t> Find(const Point2D &point) const
  {
    // Create empty list of bounding box indices
    std::vector<size_t> indices;

    // Recursively search bounding box tree for collisions
    FindRecursive(indices, *this, point, Nodes.size() - 1);

    return indices;
  }

  /// Find indices of bounding box collisions between this
  /// tree and given tree.
  ///
  /// @param tree The other tree
  /// @return Array of pairwise collisions
  std::vector<std::pair<size_t, size_t>>
  Find(const BoundingBoxTree2D &tree) const
  {
    // Create empty list of bounding box indices
    std::vector<std::pair<size_t, size_t>> indices;

    // Recursively search bounding box tree for collisions
    FindRecursive(indices, *this, tree, Nodes.size() - 1,
                  tree.Nodes.size() - 1);

    return indices;
  }

  /// Check if bounding box is empty
  bool Empty() const { return Nodes.empty(); }

  /// Clear all data
  void Clear() { Nodes.clear(); }

  /// Pretty-print
  std::string __str__() const override
  {
    return "2D bounding box tree with " + str(Nodes.size()) + " nodes";
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

  // Build bounding box tree (recursive call). The input arguments are
  // the original full vector of bounding boxes and two iterators marking
  // the beginning and end of an array of indices (into the original vector
  // of bounding boxes) to be partitioned.
  int BuildRecursive(const std::vector<BoundingBox2D> &bboxes,
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
      node.bbox = ComputeBoundingBox(bboxes, begin, end);

      // Compute main axis of bounding box
      size_t mainAxis = ComputeMainAxis(node.bbox);

      // Split boxes into two groups based on a partial sort along main axis
      auto middle = begin + (end - begin) / 2;
      if (mainAxis == 0)
        std::nth_element(begin, middle, end, LessThanX(bboxes));
      else
        std::nth_element(begin, middle, end, LessThanY(bboxes));

      // Call recursively for both groups
      node.first = BuildRecursive(bboxes, begin, middle);
      node.second = BuildRecursive(bboxes, middle, end);
    }

    // Add node to tree
    Nodes.push_back(node);

    // Return index of current node
    return Nodes.size() - 1;
  }

  // Find collisions between tree and point (recursive call)
  static void FindRecursive(std::vector<size_t> &indices,
                            const BoundingBoxTree2D &tree,
                            const Point2D &point,
                            size_t nodeIndex)
  {
    // Get current node
    const Node &node = tree.Nodes[nodeIndex];

    // Check if node and point collide (if not, do nothing)
    if (Geometry::BoundingBoxContains2D(node.bbox, point))
    {
      // Check if leaf
      if (node.IsLeaf())
      {
        // Add node index (we know that node collides with point)
        indices.push_back(node.second);
      }
      else
      {
        // Not a leaf node so check child nodes
        FindRecursive(indices, tree, point, node.first);
        FindRecursive(indices, tree, point, node.second);
      }
    }
  }

  // Find collisions between tree A and tree B (recursive call)
  static void FindRecursive(std::vector<std::pair<size_t, size_t>> &indices,
                            const BoundingBoxTree2D &treeA,
                            const BoundingBoxTree2D &treeB,
                            size_t nodeIndexA,
                            size_t nodeIndexB)
  {
    // Get current nodes
    const Node &nodeA = treeA.Nodes[nodeIndexA];
    const Node &nodeB = treeB.Nodes[nodeIndexB];

    // Check if nodes collide (if not, do nothing)
    if (Geometry::Intersect2D(nodeA.bbox, nodeB.bbox))
    {
      // Check if both nodes are leaves
      if (nodeA.IsLeaf() && nodeB.IsLeaf())
      {
        // Add node indices (we know that the nodes collide)
        indices.push_back(std::make_pair(nodeA.second, nodeB.second));
      }
      else if (nodeA.IsLeaf())
      {
        // If A is leaf, then descend B
        FindRecursive(indices, treeA, treeB, nodeIndexA, nodeB.first);
        FindRecursive(indices, treeA, treeB, nodeIndexA, nodeB.second);
      }
      else if (nodeB.IsLeaf())
      {
        // If B is leaf, then descend A
        FindRecursive(indices, treeA, treeB, nodeA.first, nodeIndexB);
        FindRecursive(indices, treeA, treeB, nodeA.second, nodeIndexB);
      }
      else
      {
        // If neither node is a leaf, traverse largest subtree
        if (nodeIndexA > nodeIndexB)
        {
          FindRecursive(indices, treeA, treeB, nodeA.first, nodeIndexB);
          FindRecursive(indices, treeA, treeB, nodeA.second, nodeIndexB);
        }
        else
        {
          FindRecursive(indices, treeA, treeB, nodeIndexA, nodeB.first);
          FindRecursive(indices, treeA, treeB, nodeIndexA, nodeB.second);
        }
      }
    }
  }

  // Compute bounding box of bounding boxes
  BoundingBox2D ComputeBoundingBox(const std::vector<BoundingBox2D> &bboxes,
                                   const std::vector<size_t>::iterator &begin,
                                   const std::vector<size_t>::iterator &end)
  {
    // Initialize bounding box
    constexpr double max = std::numeric_limits<double>::max();
    BoundingBox2D boundingBox;
    boundingBox.P = Point2D(max, max);
    boundingBox.Q = Point2D(-max, -max);

    // Iterate over bounding boxes to compute bounds
    for (auto it = begin; it != end; it++)
    {
      const BoundingBox2D &bbox = bboxes[*it];
      boundingBox.P.x = std::min(boundingBox.P.x, bbox.P.x);
      boundingBox.P.y = std::min(boundingBox.P.y, bbox.P.y);
      boundingBox.Q.x = std::max(boundingBox.Q.x, bbox.Q.x);
      boundingBox.Q.y = std::max(boundingBox.Q.y, bbox.Q.y);
    }

    return boundingBox;
  }

  // Compute main axis of bounding box
  size_t ComputeMainAxis(const BoundingBox2D &bbox)
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

}

#endif
