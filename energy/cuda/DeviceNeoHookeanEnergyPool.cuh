#pragma once

#include "common/cuda/DeviceTombstonePool.cuh"
#include "energy/elastic/NeoHookeanEnergyDefs.hpp"

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

    NeoHookeanEnergyPoolView view() const
    {
        NeoHookeanEnergyPoolView v;
        DeviceTombstonePoolView::view(v);
        v.data = data.data();
        return v;
    }

    void upload(const NeoHookeanEnergyPool& src);
};

} // namespace Energy