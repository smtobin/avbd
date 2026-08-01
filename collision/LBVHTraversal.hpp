#pragma once

#include "collision/LBVH.hpp"
#include "collision/AABB.hpp"
#include "common/WorkerThreadContext.hpp"

#include <queue>

namespace Collision
{

struct LBVHTraversal
{

    /** List of traversal tasks generated to be executed in parallel */
    std::vector<std::pair<unsigned, unsigned>> _traversal_tasks;

    /** Index of next available traversal task to be completed */
    std::atomic<unsigned> _task_counter;

    /** Offsets for merging potential collisions detected in parallel */
    std::vector<unsigned> _merge_offsets;

    /** Collides the tree with itself, starting at the designated root, and returning any overlapping leaf nodes.
     * @param bvh the LBVH to traverse
     * @param start_a the node index of one subtree root to start at
     * @param start_b the node index of the other subtree root to start at
     * @param leaf_pairs (output) the detected collisions between leaf AABBs. These are pairs of leaf indices, which correspond to the SORTED order in the collision pool.
     */
    void traverseSelfIterative(const LBVH& bvh, unsigned start_a, unsigned start_b,
                            std::vector<std::pair<unsigned,unsigned>>& leaf_pairs)
    {
        // stack of node pairs to test against each other
        // (avoids recursion)
        std::vector<std::pair<unsigned,unsigned>> work_stack;
        work_stack.reserve(256);
        work_stack.push_back({start_a, start_b});

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

    void generateTraversalTasks(const LBVH& bvh, unsigned root, unsigned target_count, std::vector<std::pair<unsigned, unsigned>>& tasks)
    {
        // helper function for a work estimate
        // use the subtree sizes multiplied together as a heuristic
        auto work_estimate = [&](unsigned a, unsigned b) -> unsigned {
            unsigned sa = (bvh.leaf_count[a] > 0) ? 1 : bvh.subtree_size[a];
            unsigned sb = (bvh.leaf_count[b] > 0) ? 1 : bvh.subtree_size[b];
            return sa * sb;
        };

        // comparator for two tasks
        using Task = std::pair<unsigned, unsigned>;
        auto cmp = [&](const Task& t1, const Task& t2) {
            return work_estimate(t1.first, t1.second) < work_estimate(t2.first, t2.second);
        };

        // max-heap, we will keep splitting the most expensive tasks until we reach our target number of tasks
        std::priority_queue<Task, std::vector<Task>, decltype(cmp)> pq(cmp);
        pq.push({root, root});

        while (pq.size() < target_count)
        {
            auto [a,b] = pq.top();
            bool leaf_a = bvh.leaf_count[a] > 0;
            bool leaf_b = bvh.leaf_count[b] > 0;

            if (leaf_a && leaf_b)
                break;  //largest remaining task can't be split further - stop

            pq.pop();

            // if the nodes are the same, split this node into 3
            if (a == b)
            {
                unsigned l = bvh.left[a];
                unsigned r = bvh.right[a];
                pq.push({l,l});
                pq.push({r,r});
                pq.push({l,r});
            }
            // if A is a leaf, or subtree B is larger than subtree A, split B
            else if (leaf_a || (!leaf_b && bvh.subtree_size[b] > bvh.subtree_size[a]))
            {
                pq.push({a, bvh.left[b]});
                pq.push({a, bvh.right[b]});
            }
            // if B is a leaf, or subtree A is larger than subtree B, split A
            else
            {
                pq.push({bvh.left[a], b});
                pq.push({bvh.right[a], b});
            }
        }

        // once we get here, we have reached our target count or all the tasks are comparing leaf nodes
        // add the tasks that remain in the priority queue to the task list
        tasks.clear();
        while (!pq.empty()) 
        { 
            tasks.push_back(pq.top()); 
            pq.pop(); 
        }
    }

    void collisionBroadPhase(WorkerThreadContext& w_ctx, const LBVH& lbvh)
    {
        // main thread generates the traversal tasks for the other threads
        if (w_ctx.idx == 0)
        {
            generateTraversalTasks(lbvh, lbvh.root, WorkerThreadContext::NUM_THREADS * 6, _traversal_tasks);
            _task_counter.store(0, std::memory_order_relaxed);
        }
        w_ctx.barrier->arrive_and_wait();

        auto& potential_collisions = w_ctx.potential_collisions;
        potential_collisions.clear();

        // fetch tasks off the task list and perform self traversal on the task subtree
        unsigned idx;
        while ((idx = _task_counter.fetch_add(1, std::memory_order_relaxed)) < _traversal_tasks.size())
        {
            auto [a,b] = _traversal_tasks[idx];
            traverseSelfIterative(lbvh, a, b, potential_collisions);
        }

        w_ctx.barrier->arrive_and_wait();
    }

    void mergeLeafPairs(WorkerThreadContext& w_ctx, const std::vector<WorkerThreadContext>& all_worker_contexts, std::vector<std::pair<unsigned, unsigned>>& merged_potential_collisions)
    {
        // main thread compute per-thread offsets
        if (w_ctx.idx == 0)
        {
            _merge_offsets.resize(WorkerThreadContext::NUM_THREADS + 1);
            _merge_offsets[0] = 0;
            for (unsigned t = 0; t < WorkerThreadContext::NUM_THREADS; t++)
            {
                _merge_offsets[t+1] = _merge_offsets[t] + all_worker_contexts[t].potential_collisions.size();
            }
            merged_potential_collisions.resize(_merge_offsets.back());
        }
        w_ctx.barrier->arrive_and_wait();

        // each thread copies its own buffer into its assigned slice
        std::copy(
            w_ctx.potential_collisions.begin(), w_ctx.potential_collisions.end(),
            merged_potential_collisions.begin() + _merge_offsets[w_ctx.idx]
        );
        w_ctx.barrier->arrive_and_wait();
    }

};

} // namespace Collision