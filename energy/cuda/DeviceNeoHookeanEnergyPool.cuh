#pragma once

#include "common/cuda/DevicePoolBase.cuh"
#include "energy/elastic/NeoHookeanEnergyDefs.hpp"
#include "energy/elastic/NeoHookeanEnergyPool.hpp"

namespace Energy
{

struct DeviceNeoHookeanEnergyPool : DevicePoolBase
{
    DeviceBuffer<NeoHookeanEnergyInfo> data;
    
    explicit DeviceNeoHookeanEnergyPool(unsigned capacity)
        : DevicePoolBase(capacity)
        , data(capacity)
    {

    }

    NeoHookeanEnergyPoolView view()
    {
        NeoHookeanEnergyPoolView v;
        DevicePoolBase::view(v);
        v.data = data.data();
        return v;
    }

    void upload(const NeoHookeanEnergyPool& src)
    {
        DevicePoolBase::upload(src);
        data.upload(src.data.data(), src.data.size());
    }
};

} // namespace Energy