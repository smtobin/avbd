#include "common/cuda/DeviceTombstonePool.cuh"
#include "common/TombstonePool.hpp"

void DeviceTombstonePool::upload(const TombstonePool& src)
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