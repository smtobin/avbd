#pragma once

#include "common/common.hpp"
#include "common/TombstonePool.hpp"
#include "common/cuda/DeviceBuffer.cuh"

/** A set of counters/flags attached to each pool.
 * One instance per pool, stored in device memory.
 */
struct PoolCounters
{
    unsigned free_count;        // number of available slots on the free stack
    unsigned high_water;        // slots [0, high_water) have been used at some point
    unsigned active_count;      // number of live elements
    unsigned error;             // sticky flags - host will check and clear at frame boundaries
};

////////////////////////////////////////////////////////////////////

/** Represents an allocation request to the TombstonePool. */
struct AllocTicket
{
    unsigned free_base;         // free_count beefore this allocation
    unsigned hw_base;           // high_water before this allocation
    unsigned from_free;         // ranks [0, from_free) are served from the free stack
    unsigned total;             // total slots granted
};

/////////////////////////////////////////////////////////////////////

/** Non-owning view into DeviceTombstonePool */
struct DeviceTombstonePoolView
{
    uint8_t* active;
    unsigned* generation;
    unsigned* free_list;

    PoolCounters* counters;
    AllocTicket* ticket;

    unsigned capacity;
};

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
    void view(DeviceTombstonePoolView& view) const
    {
        view.active = active.data();
        view.counters = counters.data();
        view.ticket = ticket.data();

        view.capacity = capacity;
    }

    void upload(const TombstonePool& src)
    {
        if (src.capacity != capacity)
            throw std::runtime_error("Host and device pool capacities do not match");

        active.upload(src.active.data(), capacity);
        // generation.upload(src.generation.data(), capacity);
        // free_list.upload(src.slots.data(), capacity);

        // counters.data()->active_count = src.count;
        // counters.data()->free_count = capacity - src.count;
        // counters.data()->high_water = src.highest_index;

        // ticket.data()->free_base = counters.data()->free_count;
        // ticket.data()->hw_base = counters.data()->high_water;
    }
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