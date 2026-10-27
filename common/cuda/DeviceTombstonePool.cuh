#pragma once

#include "common/common.hpp"
#include "common/TombstonePoolView.hpp"
#include "common/cuda/DeviceBuffer.cuh"

// forward declaration of CPU pool type
struct TombstonePool;

//////////////////////////////////////////////////////////////////

/** GPU mirror of TombstonePool - responsible for memory ownership on device side.
 */
struct DeviceTombstonePool
{
    DeviceBuffer<uint8_t> active;   // whether elements are valid

    /** Empty slot management as a stack
     * For now, don't include this in the GPU implementation
     * We will assume that all pools are static.
     * TODO: (10/16/26) figure out thread-safe way to add new elements to pool on device side
     */
    DeviceBuffer<unsigned> generation;
    DeviceBuffer<unsigned> free_list;
    DeviceBuffer<PoolCounters> counters;
    DeviceBuffer<AllocTicket> ticket;

    unsigned capacity;              // total capacity (allocated memory)

    /** Constructor initializes memory
     * @param capacity : the capacity of the memory pool
     */
    explicit DeviceTombstonePool(unsigned capacity)
        : active(capacity)
        , generation(capacity)
        , free_list(capacity)
        , counters(1)
        , ticket(1)
        , capacity(capacity)
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
    void view(TombstonePoolView& view) const
    {
        view.active = active.data();
        view.counters = counters.data();
        view.ticket = ticket.data();

        view.capacity = capacity;
    }

    void upload(const TombstonePool& src);
};

///////////////////////////////////////////////////////////////////////////////////////////////////

// __global__ void reserveSlotsKernel(DeviceTombstonePoolView v, const unsigned* d_total)
// {
//     unsigned total = *d_total;
//     PoolCounters& c = *v.counters;

//     unsigned from_free = min(total, c.free_count);
//     unsigned fresh = total - from_free;

//     if (c.error || c.high_water + fresh > v.capacity)
//     {
//         c.error |= ERR_POOL_FULL;
//         *v.ticket = {0, 0, 0, 0};
//         return;
//     }

//     *v.ticket = {c.free_count, c.high_water, from_free, total};
//     c.free_count -= from_free;
//     c.high_water += fresh;
//     c.active_count += total;
// }