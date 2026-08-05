#pragma once

#include "common/common.hpp"

#include <vector>

namespace Collision
{

/** Stores LBVH nodes in a SoA format.
 * 
 * Node index: [0, n-1] are leaves, [n, 2n-2] are internal nodes
 * 
 * In vanilla LBVH, every leaf is exactly one primitive so leaf_count is always 1. (or 0 for internal nodes)
 * Here we support leaf collapsing - merging small subtrees into a single leaf to reduce tree depth and traversal overhead.
 */
struct LBVH
{
    static constexpr unsigned INVALID = std::numeric_limits<unsigned>::max();
    
    std::vector<Real> min_x, min_y, min_z;  // bounding box min coords
    std::vector<Real> max_x, max_y, max_z;  // bounding box max coords
    std::vector<unsigned> left, right;      // left and right children
    std::vector<unsigned> parent;           // parent
    std::vector<unsigned> leaf_start;       // index into sorted order where this leaf's primitives begin
    std::vector<unsigned> leaf_count;       // how many consecutive primitives belong to this leaf
    std::vector<unsigned> subtree_size;     // total leaves under this node (leaves included)
    unsigned root;

    std::vector<std::atomic<uint8_t>> visited;  // helps track which nodes are visited during refit pass

    void resize(unsigned num_primitives)
    {
        unsigned num_nodes = 2 * num_primitives - 1;
        min_x.resize(num_nodes); min_y.resize(num_nodes); min_z.resize(num_nodes);
        max_x.resize(num_nodes); max_y.resize(num_nodes); max_z.resize(num_nodes);
        left.resize(num_nodes);  right.resize(num_nodes);
        parent.resize(num_nodes);
        leaf_start.resize(num_nodes, 0);
        leaf_count.resize(num_nodes, 0);
        subtree_size.resize(num_nodes, 0);
        if (visited.size() != num_nodes)
            visited = std::vector<std::atomic<uint8_t>>(num_nodes); // use assignment since atomics cannot be relocated if vector resizes
    }

    unsigned numPrimitives() const { return (parent.size() + 1)/2; }

    void printTree() const
    {
        _printTreeImpl(root);
    }

    void printTreeWithInfo(const CollisionPrimitivePool& pool) const
    {
        _printTreeWithInfoImpl(root, pool);
    }

    std::vector<unsigned> nodeDepths() const;

private:
    void _printTreeImpl(unsigned node, const std::string& prefix = "", bool is_left = true) const
    {
        std::string node_str;
        if (node >= numPrimitives())
            node_str = "I" + std::to_string(node - numPrimitives());
        else
            node_str = "L" + std::to_string(node);

        std::cout << prefix
                << (is_left ? "├── " : "└── ")
                << node_str << '\n';

        if (leaf_count[node] > 0) return;

        _printTreeImpl(left[node],
                prefix + (is_left ? "│   " : "    "),
                true);

        _printTreeImpl(right[node],
                prefix + (is_left ? "│   " : "    "),
                false);
    }

    void _printTreeWithInfoImpl(unsigned node, const CollisionPrimitivePool& pool, const std::string& prefix = "", bool is_left = true) const;
};

} // namespace Collision