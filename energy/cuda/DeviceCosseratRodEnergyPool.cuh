#pragma once

#include "common/cuda/DevicePoolBase.cuh"
#include "energy/elastic/CosseratRodEnergyDefs.hpp"
#include "energy/elastic/CosseratRodEnergyPool.hpp"

namespace Energy
{

struct DeviceCosseratRodEnergyPool : DevicePoolBase
{
    DeviceBuffer<CosseratRodEnergyInfo> data;
    
    explicit DeviceCosseratRodEnergyPool(unsigned capacity)
        : DevicePoolBase(capacity)
        , data(capacity)
    {

    }

    CosseratRodEnergyPoolView view()
    {
        CosseratRodEnergyPoolView v;
        DevicePoolBase::view(v);
        v.data = data.data();

        return v;
    }

    void upload(const CosseratRodEnergyPool& src)
    {
        DevicePoolBase::upload(src);
        data.upload(src.data.data(), src.data.size());
    }
};

} // namespace Energy