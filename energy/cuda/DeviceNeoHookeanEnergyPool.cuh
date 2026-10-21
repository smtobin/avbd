#pragma once

#include "common/cuda/Mat3.cuh"
#include "common/cuda/DeviceTombstonePool.cuh"
#include "energy/elastic/NeoHookeanEnergyPool.hpp"

namespace Energy
{


struct DeviceNeoHookeanEnergyPoolView : DeviceTombstonePoolView
{
    NeoHookeanEnergyInfo* data;
};

struct DeviceNeoHookeanEnergyPool : DeviceTombstonePool
{
    DeviceBuffer<NeoHookeanEnergyInfo> data;
    
    explicit DeviceNeoHookeanEnergyPool(unsigned capacity)
        : DeviceTombstonePool(capacity)
        , data(capacity)
    {

    }

    DeviceNeoHookeanEnergyPoolView view() const
    {
        DeviceNeoHookeanEnergyPoolView v;
        DeviceTombstonePoolView::view(v);
        v.data = data.data();

        return v;
    }

    void upload(const NeoHookeanEnergyPool& src)
    {
        DeviceTombstonePool::upload(src);
        data.upload(src.data);
    }
};

} // namespace Energy