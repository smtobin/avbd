#pragma once

#include "common/common.hpp"
#include "common/PoolBase.hpp"
#include "common/cuda/DeviceValue.cuh"
#include "common/cuda/DeviceBuffer.cuh"

/** Append-only base for all GPU memory pools in the sim. */
struct DevicePoolBase
{
    unsigned capacity;
    unsigned host_high_water;

    DeviceValue<unsigned> high_water;

    DevicePoolBase(unsigned capacity_)
        : capacity(capacity_)
        , host_high_water(0)
        , high_water(0)
    {}

    void view(PoolBaseView& v)
    {
        v.capacity = capacity;
        v.high_water = high_water.data();
    }

    void upload(const PoolBase& src)
    {
        if (src.capacity != capacity)
            throw std::runtime_error("Host and device pool capacities do not match");

        host_high_water = src.high_water;
        high_water.set(src.high_water);
    }

    void download(PoolBase& dst) const
    {
        if (dst.capacity != capacity)
            throw std::runtime_error("Host and device pool capacities do not match");

        dst.high_water = high_water.get();
    }
};