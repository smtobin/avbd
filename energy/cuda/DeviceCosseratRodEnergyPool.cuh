#pragma once

#include "common/cuda/DeviceTombstonePool.cuh"
#include "energy/elastic/CosseratRodEnergyDefs.hpp"
#include "energy/elastic/CosseratRodEnergyPool.hpp"

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

    CosseratRodEnergyPoolView view()
    {
        CosseratRodEnergyPoolView v;
        DeviceTombstonePool::view(v);
        v.data = data.data();

        return v;
    }

    void upload(const CosseratRodEnergyPool& src)
    {
        DeviceTombstonePool::upload(src);
        data.upload(src.data.data(), src.data.size());
    }
};

} // namespace Energy