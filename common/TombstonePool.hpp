#pragma once

#include "common/common.hpp"
#include "common/TombstonePoolView.hpp"

#include <vector>
#include <atomic>
#include <numeric>

/** Base struct for dynamic structure-of-array memory "pools" which have elements added and removed throughout the sim.
 * 
 * Rather than constantly shifting elements to eliminate gaps, keeps a list of "tombstones" so that we know where the gaps are.
 */
struct TombstonePool
{
    
    std::vector<uint8_t> active;   // whether elements are valid
    unsigned capacity;          // total capacity (allocated memory)
    // unsigned count;             // current number of active elements
    // unsigned highest_index;     // highest active index

    PoolCounters counters;

    /** Empty slot management as a stack */
    std::vector<unsigned> slots;    // stack of empty slots in the pool
    std::atomic<int> top;   // points to the "top" of the stack - i.e. the first valid (free) slot in the pool

    /** Store a generation per slot. Helps us know what slots have been updated. */
    std::vector<unsigned> generation;

    /** Constructor initializes memory
     * @param capacity : the capacity of the memory pool
     */
    explicit TombstonePool(unsigned capacity)
        : active(capacity, false)
        , capacity(capacity)
        // , count(0)
        // , highest_index(0)
        , counters({capacity, 0, 0, 0, -1})
        , slots(capacity)
        , top(0)
        , generation(capacity)
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

    /** Returns a view of the TombstonePool
     * TODO: (10/27/26) This is more to align CPU with GPU - when would we use TombstonePoolView on CPU?
     * For now, PoolCounters and AllocTicket are not set in the view.
     */
    void view(TombstonePoolView& v)
    {
        v.active = active.data();
        v.generation = generation.data();
        v.free_list = slots.data();
        v.counters = &counters;

        v.capacity = capacity;
    }

    /** The total size of the pool, from index 0 to the highest active index */
    int totalSize() const { return counters.highest_index+1; }

    /** Allocates space for a new element based on the free slots.
     * @returns the index of the slot that has been allocated for the new element.
     */
    unsigned allocSlot()
    {
        int idx = top.fetch_add(1);
        if (idx >= (int)capacity)
        {
            top.fetch_sub(1);   // undo add
            throw std::runtime_error("Pool has exceeded capacity");
        }

        unsigned slot = slots[idx];
        active[slot] = true;
        generation[slot]++;
        counters.active_count++;
        counters.highest_index = std::max(counters.highest_index, static_cast<int>(slot));

        // std::cout << "Alloc slot..." << std::endl;
        // std::cout << "  New highest index: " << highest_index << std::endl;
        // std::cout << "  Slot: " << slot << std::endl;

        return slot;
    }

    /** Frees space for removing an element from a slot.
     * @param slot : the position of the slot to free
     */
    void freeSlot(unsigned slot)
    {
        // slot must be active to free it
        if (!active[slot])
        {
            throw std::runtime_error("freeSlot called on inactive slot");
        }

        // update highest index if the slot being freed is the highest index
        if (static_cast<int>(slot) == counters.highest_index)
        {
            while (counters.highest_index >= 0 &&
                !active[counters.highest_index])
            {
                --counters.highest_index;
            }
        }

        active[slot] = false;
        counters.active_count--;

        int idx = top.fetch_sub(1) - 1;
        slots[idx] = slot;
    }

    /** Iterator subclass that iterates through active indices. */
    class iterator
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = unsigned;
        using difference_type = std::ptrdiff_t;
        using pointer = const unsigned*;
        using reference = unsigned;

        iterator(const TombstonePool* pool, unsigned index)
            : _pool(pool)
            , _index(index)
        {
            skipInactive();
        }

        unsigned operator*() const
        {
            return _index;
        }

        iterator& operator++()
        {
            ++_index;
            skipInactive();
            return *this;
        }

        iterator operator++(int)
        {
            iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const iterator& other) const
        {
            return _pool == other._pool &&
                   _index == other._index;
        }

        bool operator!=(const iterator& other) const
        {
            return !(*this == other);
        }

    private:
        void skipInactive()
        {
            const unsigned limit = _pool->counters.highest_index + 1;

            while (_index < limit &&
                   !_pool->active[_index])
            {
                ++_index;
            }
        }

        const TombstonePool* _pool;
        unsigned _index;
    };

    iterator begin() const
    {
        return iterator(this, 0);
    }

    iterator end() const
    {
        return iterator(this, counters.highest_index + 1);
    }
};

