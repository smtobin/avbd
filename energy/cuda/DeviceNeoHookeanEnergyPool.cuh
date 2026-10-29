#pragma once

#include "common/cuda/DeviceTombstonePool.cuh"
#include "energy/elastic/NeoHookeanEnergyDefs.hpp"
#include "energy/elastic/NeoHookeanEnergyPool.hpp"

namespace Energy
{

struct DeviceNeoHookeanEnergyPool : DeviceTombstonePool
{
    DeviceBuffer<NeoHookeanEnergyInfo> data;
    
    explicit DeviceNeoHookeanEnergyPool(unsigned capacity)
        : DeviceTombstonePool(capacity)
        , data(capacity)
    {

    }

    NeoHookeanEnergyPoolView view()
    {
        NeoHookeanEnergyPoolView v;
        DeviceTombstonePool::view(v);
        v.data = data.data();
        return v;
    }

    void upload(const NeoHookeanEnergyPool& src)
    {
        DeviceTombstonePool::upload(src);
        data.upload(src.data.data(), src.data.size());
    }
};

} // namespace Energy