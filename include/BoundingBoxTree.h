// Axis-aligned bounding box tree (AABB tree)
// Copyright (C) 2020 Anders Logg

#ifndef DTCC_BOUNDING_BOX_TREE_H
#define DTCC_BOUNDING_BOX_TREE_H

#include <vector>
#include <limits>

#include "Point.h"

namespace DTCC
{

  // Note: This is a new implementation of the corresponding algorithm
  // in FEniCS/DOLFIN (by the same author) using simpler data structures
  // and simpler design (separate 2D and 3D implementations).

  // Axis-aligned bouding box tree (2D).
  class BoundingBoxTree2D
  {
  public:

    // Tree node data structure. Each has two children (unless it is a
    // leaf node) and a bounding box defining the region of the node
    // A leaf-node is indicated by setting the first child node to -1.
    struct Node
    {
      // Index to first child node
      int first{};

      // Index to second child node
      int second{};

      // Bounding box of node
      BoundingBox2D bbox;
    };

    // Nodes of the bounding box tree. The tree is built bottom-up, meaning
    // that the leaf nodes will be added first and the top node last.
    std::vector<Node> Nodes;

    // Build bounding box tree for bounding boxes
    void Build(const std::vector<BoundingBox2D>& bboxes)
    {
      std::cout << "BoundingTree: Building 2D bounding box tree..." << std::endl;

      // Initialize indices of bounding boxes to be sorted
      std::vector<size_t> indices(bboxes.size());
      for (size_t i = 0; i < bboxes.size(); i++)
        indices[i] = i;

      // Clear tree if built before
      Nodes.clear();

      // Recursively build bounding box tree
      BuildRecursive(bboxes, indices.begin(), indices.end());

      //
    }

    // Find indices of bounding boxes containing point
    std::vector<size_t> Find(const Point2D& point) const
    {
      // Create empty list of bounding box indices
      std::vector<size_t> indices;

      // Recursively search bounding box tree for collisions
      FindRecursive(indices, point, Nodes.size() - 1);

      return indices;
    }

  private:

    /// Comparison for partial sorting of bounding boxes along x-axis
    struct LessThanX
    {
      const std::vector<BoundingBox2D>& bboxes;
      LessThanX(const std::vector<BoundingBox2D>& bboxes): bboxes(bboxes) {}
      bool operator()(size_t i, size_t j)
      {
        return bboxes[i].P.x + bboxes[i].Q.x < bboxes[j].P.x + bboxes[j].Q.x;
      }
    };

    /// Comparison for partial sorting of bounding boxes along y-axis
    struct LessThanY
    {
      const std::vector<BoundingBox2D>& bboxes;
      LessThanY(const std::vector<BoundingBox2D>& bboxes): bboxes(bboxes) {}
      bool operator()(size_t i, size_t j)
      {
        return bboxes[i].P.y + bboxes[i].Q.y < bboxes[j].P.y + bboxes[j].Q.y;
      }
    };

    // Build bounding box tree (recursive call). The input arguments are
    // the original full vector of bounding boxes and two iterators marking
    // the beginning and end of an array of indices (into the original vector
    // of bounding boxes) to be partitioned.
    int BuildRecursive(const std::vector<BoundingBox2D>& bboxes,
                       const std::vector<size_t>::iterator& begin,
                       const std::vector<size_t>::iterator& end)
    {
      // Create empty node
      Node node;

      // Check if we reached a leaf
      if (end - begin == 1)
      {
        node.bbox = bboxes[*begin];
        node.first = -1;
        node.second = -1;
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

    // Findcollisions (recursive call)
    void FindRecursive(std::vector<size_t>& indices,
                       const Point2D& point,
                       size_t nodeIndex) const
    {
      // Get current node
      const Node& node = Nodes[nodeIndex];

      // Check if node contains point (do nothing if not contained)
      if (Geometry::BoundingBoxContains2D(node.bbox, point))
      {
        if (node.first == -1)
        {
          // Leaf node containing point so add it
          indices.push_back(nodeIndex);
        }
        else
        {
          // Not a leaf node so check child nodes
          FindRecursive(indices, point, node.first);
          FindRecursive(indices, point, node.second);
        }
      }
    }

    // Compute bounding box of bounding boxes
    BoundingBox2D ComputeBoundingBox(const std::vector<BoundingBox2D>& bboxes,
                                     const std::vector<size_t>::iterator& begin,
                                     const std::vector<size_t>::iterator& end)
    {
      // Initialize bounding box
      constexpr double max = std::numeric_limits<double>::max();
      BoundingBox2D boundingBox;
      boundingBox.P = Point2D(max, max);
      boundingBox.Q = Point2D(-max, -max);

      // Iterate over bounding boxes to compute bounds
      for (auto it = begin; it != end; it++)
      {
        const BoundingBox2D& bbox = bboxes[*it];
        boundingBox.P.x = std::min(boundingBox.P.x, bbox.P.x);
        boundingBox.P.y = std::min(boundingBox.P.y, bbox.P.y);
        boundingBox.Q.x = std::max(boundingBox.Q.x, bbox.Q.x);
        boundingBox.Q.y = std::max(boundingBox.Q.y, bbox.Q.y);
      }

      return boundingBox;
    }

    // Compute main axis of bounding boxes
    size_t ComputeMainAxis(const BoundingBox2D& bbox)
    {
      const double dx = bbox.Q.x - bbox.P.x;
      const double dy = bbox.Q.y - bbox.P.y;
      return (dx > dy ? 0 : 1);
    }

  };

  class BoundingBoxTree3D
  {
    // FIXME: Not yet implemented
  };

  std::ostream &operator<<(std::ostream &stream, const BoundingBoxTree2D& bbtree)
  {
    stream << "2D bounding box tree with " << bbtree.Nodes.size() << " nodes" << std::endl;
    return stream;
  }

}

#endif
