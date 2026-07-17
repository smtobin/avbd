#pragma once

#include "collision/LBVH.hpp"
#include "collision/AABB.hpp"

namespace Collision
{

struct LBVHTraversal
{

    static void traversePair(const LBVH& bvh, unsigned node_a, unsigned node_b, std::vector<std::pair<unsigned, unsigned>>& out_pairs)
    {
        AABB box_a {
            {bvh.min_x[node_a], bvh.min_y[node_a], bvh.min_z[node_a]},
            {bvh.max_x[node_a], bvh.max_y[node_a], bvh.max_z[node_a]}
        };
        AABB box_b {
            {bvh.min_x[node_b], bvh.min_y[node_b], bvh.min_z[node_b]},
            {bvh.max_x[node_b], bvh.max_y[node_b], bvh.max_z[node_b]}
        };

        // make sure boxes overlap
        if (!box_a.overlaps(box_b))
            return;

        bool leaf_a = bvh.leaf_count[node_a] > 0;
        bool leaf_b = bvh.leaf_count[node_b] > 0;

        // if both are leaves - add this pair as output for narrow-phase
        if (leaf_a && leaf_b)
        {
            if (node_a < node_b)    // avoid adding both (a,b) and (b,a), and avoid (a,a)
                out_pairs.push_back({node_a, node_b});
            return;
        }

        // traverse pair of subtrees
        // descend down whichever tree is larger
        if (leaf_a || (!leaf_b && bvh.subtree_size[node_b] > bvh.subtree_size[node_a]))
        {
            traversePair(bvh, node_a, bvh.left[node_b], out_pairs);
            traversePair(bvh, node_a, bvh.right[node_b], out_pairs);
        }
        else
        {
            traversePair(bvh, bvh.left[node_a], node_b, out_pairs);
            traversePair(bvh, bvh.right[node_a], node_b, out_pairs);
        }
    }

    static void traverseSelf(const LBVH& bvh, unsigned node, std::vector<std::pair<unsigned, unsigned>>& out_pairs)
    {
        // single leaf can't self-collide
        if (bvh.leaf_count[node] > 0)
            return;

        unsigned l = bvh.left[node];
        unsigned r = bvh.right[node];
        // left internal self-collisions
        traverseSelf(bvh, l, out_pairs);
        // right internal self-collisions
        traverseSelf(bvh, r, out_pairs);
        // cross term, left vs. right
        traversePair(bvh, l, r, out_pairs);
    }

    void traverseSelfIterative(const LBVH& bvh, unsigned root,
                            std::vector<std::pair<unsigned,unsigned>>& out_pairs)
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
                    out_pairs.push_back({node_a, node_b});
                return;
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