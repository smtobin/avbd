#pragma once

#include "common/common.hpp"
#include "common/SpinBarrier.hpp"
#include "common/TombstonePool.hpp"
#include "collision/AABB.hpp"
#include "collision/DetectedCollision.hpp"

/** Basically a memory pool for the different stages of the time step for each worker.
 * 
 * Align on 64-byte boundaries so that each WorkerThreadContext is on its own cache line.
 */
struct alignas(64) WorkerThreadContext
{
    inline static unsigned NUM_THREADS = 1;

    // thread index
    unsigned idx;
    // pointer to global spin barrier
    SpinBarrier* barrier;

    // stores potentially colliding leaf pairs detected during broad phase
    std::vector<std::pair<unsigned, unsigned>> potential_collisions; 

    // stores detected collisions from narrow phase
    std::vector<DetectedCollision> detected_collisions;

    // scratch memory for computation
    union 
    {
        // for building the BVH
        struct {
            Collision::AABB scene_box;  // this thread's AABB around all primitives
        }   BuildBVHContext;
        
        // for parallel radix sort
        struct {
            size_t local_count[256];    // this thread's histogram for the current pass
        } RadixSortContext;
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
    std::pair<unsigned, unsigned> computeStartEnd(const TombstonePool& pool) const
    {
        return computeStartEnd(pool.totalSize());
    }
};