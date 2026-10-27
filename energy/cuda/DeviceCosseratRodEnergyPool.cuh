#pragma once

#include "common/cuda/DeviceTombstonePool.cuh"
#include "energy/elastic/CosseratRodEnergyDefs.hpp"

namespace Energy
{

struct DeviceCosseratRodEnergyPool : DeviceTombstonePool
{
    DeviceBuffer<CosseratRodEnergyInfo> data;
    
    explicit DeviceCosseratRodEnergyPool(unsigned capacity)
        : DeviceTombstonePool(capacity)
        , data(capacity)
    {

    }

    CosseratRodEnergyPoolView view() const
    {
        CosseratRodEnergyPoolView v;
        DeviceTombstonePool::view(v);
        v.data = data.data();

        return v;
    }

    void upload(const CosseratRodEnergyPool& src);
};

} // namespace Energy