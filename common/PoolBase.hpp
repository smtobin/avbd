#pragma once

#include "common/common.hpp"

#include <vector>
#include <atomic>

/** Non-owning view into PoolBase. */
struct PoolBaseView
{
    unsigned* high_water;

    bool isActive(unsigned i) const { return i < *high_water; }
};

////////////////////////////////////////////////////////////////////////////////////////////////////

/** Append-only base for all memory pools in the sim.
 * 
 * Assumes that concurrent reads and appends do not happen.
 * i.e., appends happen in one phase (potentially concurrently) and reads happen in another phase.
 */
struct PoolBase
{
    /** Sentinel value indicating an invalid index. */
    static constexpr unsigned INVALID = ~0u;

    /** Capacity of the pool. Should not change during the sim. */
    unsigned capacity;

    /** High water mark. Indices [0, high_water) have been used at some point during the sim, and
     * indices [high_water, capacity) are unused.
     * 
     * Essentially provides an upper bound for the pool.
     */
    unsigned high_water;

    /** Error flag if the pool overflows. */
    uint8_t  overflowed = 0;

    explicit PoolBase(unsigned capacity_)
        : capacity(capacity_)
        , high_water(0)
    {}

    /** Returns a non-owning view into the pool */
    void view(PoolBaseView& v)
    {
        v.high_water = &high_water;
    }

    /** Reserve n contiguous slots. Returns the base index, or INVALID if full.
     * Assumes that reads on newly allocated slots are not happening concurrently.
     * (otherwise it may read data that is currently being copied)
     */
    unsigned allocate(unsigned n = 1)
    {
        std::atomic_ref<unsigned> hw(high_water);
        unsigned cur = hw.load(std::memory_order_relaxed);
        do {
            if (n > capacity - cur) {            // cur <= capacity always holds, so no underflow
                std::atomic_ref<uint8_t>(overflowed).store(1, std::memory_order_relaxed);
                return INVALID;
            }
        } while (!hw.compare_exchange_weak(cur, cur + n,
                                           std::memory_order_relaxed,
                                           std::memory_order_relaxed));
        return cur;
    }

    /** Uniform interface so generic loops work on either pool kind. */
    bool isActive(unsigned i) const { return i < high_water; }

    /** NOT thread-safe: call between phases from a single thread. */
    void reset()
    {
        high_water = 0;
        overflowed = 0;
    }

protected:
    ~PoolBase() = default;   // not meant to be deleted through a base pointer
};