#pragma once

#include "collision/LBVH.hpp"
#include "collision/AABB.hpp"

namespace Collision
{

struct LBVHTraversal
{
    /** Collides the tree with itself, starting at the designated root, and returning any overlapping leaf nodes.
     * @param bvh the LBVH to traverse
     * @param root the node index in the LBVH to start at
     * @param leaf_pairs (output) the detected collisions between leaf AABBs. These are pairs of leaf indices, which correspond to the SORTED order in the collision pool.
     */
    static void traverseSelfIterative(const LBVH& bvh, unsigned root,
                            std::vector<std::pair<unsigned,unsigned>>& leaf_pairs)
    {
        // stack of node pairs to test against each other
        // (avoids recursion)
        std::vector<std::pair<unsigned,unsigned>> work_stack;
        work_stack.reserve(256);
        work_stack.push_back({root, root});

        // process stack until empty
        while (!work_stack.empty())
        {
            auto [node_a, node_b] = work_stack.back();
            work_stack.pop_back();

            AABB box_a{ {bvh.min_x[node_a], bvh.min_y[node_a], bvh.min_z[node_a]},
                        {bvh.max_x[node_a], bvh.max_y[node_a], bvh.max_z[node_a]} };
            AABB box_b{ {bvh.min_x[node_b], bvh.min_y[node_b], bvh.min_z[node_b]},
                        {bvh.max_x[node_b], bvh.max_y[node_b], bvh.max_z[node_b]} };

            if (!box_a.overlaps(box_b))
                continue;

            bool leaf_a = bvh.leaf_count[node_a] > 0;
            bool leaf_b = bvh.leaf_count[node_b] > 0;

            // if both are leaves - add this pair as output for narrow-phase
            if (leaf_a && leaf_b)
            {
                if (node_a < node_b)    // avoid adding both (a,b) and (b,a), and avoid (a,a)
                    leaf_pairs.push_back({node_a, node_b});
                continue;
            }

            // if both nodes are the same, we are colliding this subtree against itself
            // add 3 checks: left-left, right-right, right-left
            if (node_a == node_b)
            {
                unsigned l = bvh.left[node_a], r = bvh.right[node_a];
                work_stack.push_back({l, l});
                work_stack.push_back({r, r});
                work_stack.push_back({l, r});
                continue;
            }

            // traverse pair of subtrees
            // descend down whichever tree is larger (heuristic for balance)
            if (leaf_a || (!leaf_b && bvh.subtree_size[node_b] > bvh.subtree_size[node_a]))
            {
                work_stack.push_back({node_a, bvh.left[node_b]});
                work_stack.push_back({node_a, bvh.right[node_b]});
            }
            else
            {
                work_stack.push_back({bvh.left[node_a], node_b});
                work_stack.push_back({bvh.right[node_a], node_b});
            }
        }
    }

};

} // namespace Collision