#pragma once

#include "common/common.hpp"

/** A set of counters/flags attached to each pool.
 * One instance per pool, stored in device memory.
 */
struct PoolCounters
{
    unsigned free_count;        // number of available slots on the free stack
    unsigned high_water;        // slots [0, high_water) have been used at some point
    unsigned active_count;      // number of live elements
    unsigned error;             // sticky flags - host will check and clear at frame boundaries

    int highest_index;          // adding this for synergy with CPU TombstonePool
    /** TODO: (10/27/26) Reconcile this difference and align CPU and GPU versions of TombstonePool. */
};

////////////////////////////////////////////////////////////////////

/** Represents an allocation request to the TombstonePool.
 * This is necessary for batched allocations on the GPU.
 */
struct AllocTicket
{
    unsigned free_base;         // free_count beefore this allocation
    unsigned hw_base;           // high_water before this allocation
    unsigned from_free;         // ranks [0, from_free) are served from the free stack
    unsigned total;             // total slots granted
};

/////////////////////////////////////////////////////////////////////

/** TODO: (10/27/26) There is some divergence in the behavior of CPU TombstonePool and GPU TombstonePool.
 * E.g. high_water vs. highest_index, no AllocTicket on CPU, etc.
 * Need to reconcile this (or not).
 * TombstonePoolView more crucial for GPU than for CPU.
 */
/** Non-owning view into TombstonePool */
struct TombstonePoolView
{
    uint8_t* active;
    unsigned* generation;
    unsigned* free_list;

    PoolCounters* counters;
    AllocTicket* ticket;

    unsigned capacity;
};