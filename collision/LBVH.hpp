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
    std::vector<Real> min_x, min_y, min_z;  // bounding box min coords
    std::vector<Real> max_x, max_y, max_z;  // bounding box max coords
    std::vector<unsigned> left, right;      // left and right children
    std::vector<unsigned> parent;           // parent
    std::vector<unsigned> leaf_start;       // index into sorted order where this leaf's primitives begin
    std::vector<unsigned> leaf_count;       // how many consecutive primitives belong to this leaf

    void resize(unsigned num_primitives)
    {
        unsigned num_nodes = 2 * num_primitives - 1;
        min_x.resize(num_nodes); min_y.resize(num_nodes); min_z.resize(num_nodes);
        max_x.resize(num_nodes); max_y.resize(num_nodes); max_z.resize(num_nodes);
        left.resize(num_nodes);  right.resize(num_nodes);
        parent.resize(num_nodes);
        leaf_start.resize(num_nodes, 0);
        leaf_count.resize(num_nodes, 0);
    }
};

} // namespace Collision