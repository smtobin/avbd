#pragma once

#include "common/cuda/DeviceTombstonePool.cuh"
#include "energy/elastic/CosseratRodEnergyPool.hpp"

namespace Energy
{

struct DeviceCosseratRodEnergyPoolView : DeviceTombstonePoolView
{
    CosseratRodEnergyInfo* data;
};

struct DeviceCosseratRodEnergyPool
{
    DeviceBuffer<CosseratRodEnergyInfo> data;
    
    explicit DeviceCosseratRodEnergyPool(unsigned capacity)
        : DeviceTombstonePool(capacity)
        , data(capacity)
    {

    }

    DeviceCosseratRodEnergyPoolView view() const
    {
        DeviceCosseratRodEnergyPoolView v;
        DeviceTombstonePoolView::view(v);
        v.data = data.data();

        return v;
    }

    void upload(const CosseratRodEnergyPool& src)
    {
        DeviceTombstonePool::upload(src);
        data.upload(src.data);
    }
};

} // namespace Energy