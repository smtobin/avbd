#pragma once

#include "common/common.hpp"

#include <vector>
#include <atomic>

/** Base struct for dynamic structure-of-array memory "pools" which have elements added and removed throughout the sim.
 * 
 * Rather than constantly shifting elements to eliminate gaps, keeps a list of "tombstones" so that we know where the gaps are.
 */
struct TombstonePool
{
    
    std::vector<bool> active;   // whether elements are valid
    unsigned capacity;          // total capacity (allocated memory)
    unsigned count;             // current number of active elements
    unsigned highest_index;     // highest active index

    /** Empty slot management as a stack */
    std::vector<unsigned> slots;    // stack of empty slots in the pool
    std::atomic<int> top;   // points to the "top" of the stack - i.e. the first valid (free) slot in the pool

    /** Constructor initializes memory
     * @param capacity : the capacity of the memory pool
     */
    explicit TombstonePool(unsigned capacity)
        : active(capacity, false)
        , capacity(capacity)
        , count(0)
        , highest_index(0)
        , slots(capacity)
        , top(0)
    {
        // free slots count up from 0 to capacity
        std::iota(slots.begin(), slots.end(), 0u);
    }

    // Non-copyable
    TombstonePool(const TombstonePool&) = delete;
    TombstonePool& operator=(const TombstonePool&) = delete;

    // Movable
    TombstonePool(TombstonePool&&) = default;
    TombstonePool& operator=(TombstonePool&&) = default;

    /** Allocates space for a new element based on the free slots.
     * @returns the index of the slot that has been allocated for the new element.
     */
    unsigned allocSlot()
    {
        int idx = top.fetch_add(1);
        if (idx >= capacity)
        {
            top.fetch_sub(1);   // undo add
            throw std::runtime_error("Pool has exceeded capacity");
        }

        unsigned slot = slots[idx];
        active[slot] = true;
        count++;
        highest_index = std::max(highest_index, slot);
        return slot;
        
        return pool.slots[idx];
    }

    /** Frees space for removing an element from a slot.
     * @param slot : the position of the slot to free
     */
    void freeSlot(unsigned slot)
    {
        if (!active[slot])
        {
            throw std::runtime_error("freeSlot called on inactive slot");
        }

        active[slot] = false;
        count--;

        int idx = top.fetch_sub(1) - 1;
        slots[idx] = slot;
    }
};

