#pragma once

#include "common/cuda/DeviceBuffer.hpp"

/** Non-owning view into DeviceTombstonePool */
struct DeviceTombstonePoolView
{
    uint8_t* active;
    unsigned* capacity;
    unsigned* count;
    unsigned* highest_index;
};

//////////////////////////////////////////////////////////////////

/** GPU mirror of TombstonePool - responsible for memory ownership on device side.
 */
struct DeviceTombstonePool
{
    DeviceBuffer<uint8_t> active;   // whether elements are valid
    unsigned capacity;              // total capacity (allocated memory)
    
    /** TODO: (10/16/26) Eventually these will have to be cuda-allocated scalars. */
    unsigned count;                 // current number of active elements
    unsigned highest_index;         // highest active index

    /** Empty slot management as a stack
     * For now, don't include this in the GPU implementation
     * We will assume that all pools are static.
     * TODO: (10/16/26) figure out thread-safe way to add new elements to pool on device side
     */
    // DeviceBuffer<unsigned> slots;    // stack of empty slots in the pool
    // std::atomic<int> top;   // points to the "top" of the stack - i.e. the first valid (free) slot in the pool

    /** Store a generation per slot. Helps us know what slots have been updated. */
    // std::vector<unsigned> generation;

    /** Constructor initializes memory
     * @param capacity : the capacity of the memory pool
     */
    explicit DeviceTombstonePool(unsigned capacity)
        : active(capacity)
        , capacity(capacity)
        , count(0)
        , highest_index(0)
    {
        // nothing here
    }

    // Non-copyable
    DeviceTombstonePool(const DeviceTombstonePool&) = delete;
    DeviceTombstonePool& operator=(const DeviceTombstonePool&) = delete;

    // Movable
    DeviceTombstonePool(DeviceTombstonePool&&) = default;
    DeviceTombstonePool& operator=(DeviceTombstonePool&&) = default;

    /** Fill out device view */
    void view(DeviceTombstonePoolView& view) const
    {
        view.active = active.data();
        view.capacity = capacity;

        /** TODO: (10/16/26) these do not change because we are not allocating new elements. So we can assign them by value. */
        view.count = count;
        view.highest_index = highest_index;
    }

    void upload(const TombstonePool& src)
    {
        if (src.capacity != capacity)
        {
            throw std::runtime_error(
                "Host and device pool capacities do not match");
        }

        active.upload(src.active.data(), capacity);

        count = src.count;
        highest_index = src.highest_index;
    }

};