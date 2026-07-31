#pragma once

#include "common/common.hpp"
#include "common/SpinBarrier.hpp"
#include "collision/AABB.hpp"

/** Basically a memory pool for the different stages of the time step for each worker.
 * 
 * Align on 64-byte boundaries so that each WorkerThreadContext is on its own cache line.
 */
struct alignas(64) WorkerThreadContext
{
    static unsigned NUM_THREADS = 1;

    // thread index
    unsigned idx;
    // pointer to global spin barrier
    SpinBarrier* barrier;

    // scratch memory for computation
    union 
    {
        // for building the BVH
        struct {
            Collision::AABB scene_box;
        }   BuildBVHContext; 
    };

    /** Constructor for easily setting up context */
    WorkerThreadContext(unsigned idx_, SpinBarrier* barrier_)
        : idx(idx_), barrier(barrier_)
    {}

    /** Compute start and end indices for worker thread when iterating over a number of objects. */
    std::pair<unsigned, unsigned> computeStartEnd(unsigned total_num) const
    {
        unsigned start = total_num * idx / NUM_THREADS;
        unsigned end   = total_num * (idx + 1) / NUM_THREADS;
        return std::make_pair(start, end);
    }

    /** Compute start and end indices for worker thread when iterating over a TombstonePool. */
    std::pair<unsigned, unsigned> computeStartEnd(unsigned w_idx, const TombstonePool& pool) const
    {
        unsigned num = pool.totalSize();
        return computeStartEnd(w_idx, pool.totalSize());
    }
};